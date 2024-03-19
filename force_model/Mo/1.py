from typing import Optional
import argparse
import copy
import csv
from scipy import stats
import numpy as np
import shutil
import time
import torch
import torch.nn as nn
import torch.optim as optim
#import wandb
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler, PowerTransformer
from torch.utils.data import Dataset, DataLoader
from torch.utils.data.dataset import random_split

from utils import AverageMeter, Normalizer, ProgressMeter


class RBFExpansion(nn.Module):
    """Expand interatomic distances with radial basis functions."""

    def __init__(
            self,
            vmin: float = 0,
            vmax: float = 8,
            bins: int = 40,
            lengthscale: Optional[float] = None,
    ):
        """Register torch parameters for RBF expansion."""
        super().__init__()
        self.vmin = vmin
        self.vmax = vmax
        self.bins = bins
        self.register_buffer(
            "centers", torch.linspace(self.vmin, self.vmax, self.bins)
        )

        if lengthscale is None:
            # SchNet-style
            # set lengthscales relative to granularity of RBF expansion
            self.lengthscale = np.diff(self.centers).mean()
            self.gamma = 1 / self.lengthscale

        else:
            self.lengthscale = lengthscale
            self.gamma = 1 / (lengthscale ** 2)

    def forward(self, distance: torch.Tensor) -> torch.Tensor:
        """Apply RBF expansion to interatomic distance tensor."""
        return torch.exp(
            -self.gamma * (distance.unsqueeze(1) - self.centers) ** 2
        )


class NNDataset(Dataset):
    def __init__(self):
        super(NNDataset, self).__init__()
        with open('all_input', 'r') as f:
            self.data = np.fromstring(f.read(), sep='\n').reshape(-1,244)



 #       self.bonds_trf = PowerTransformer()
 #       bonds = self.bonds_trf.fit_transform(bonds)


        # self.input_scaler = MinMaxScaler()
        # self.data = self.input_scaler.fit_transform(data)

        with open('all_output', 'r') as f:
            self.target = np.fromstring(f.read(), sep='\n').reshape(-1, 3)

        # self.target_scaler = MinMaxScaler()
        # self.target = self.target_scaler.fit_transform(energy)

    def __len__(self):
        return len(self.target)

    def __getitem__(self, index):
        data = self.data[index]
        target = self.target[index]

        data = torch.tensor(data, dtype=torch.float)
        target = torch.tensor(target, dtype=torch.float)
        return data, target


class NNModel(nn.Module):
    def __init__(self, input_dim):
        super(NNModel, self).__init__()

    #    self.fc1 = nn.Linear(input_dim,32)
    #    self.fc2 = nn.Linear(32, 16)
    #    self.fc3 = nn.Linear(16, 6)
        self.mc1 = nn.Linear(244,64)
        self.mc2 = nn.Linear(64,32)
        self.mc3 = nn.Linear(32,3)

    #    self.bn1 = nn.BatchNorm1d(input_dim)
        self.bn1 = nn.BatchNorm1d(244)
        self.activation = nn.ReLU()


    def forward(self, x):
        x=self.bn1(x)
        x=self.mc1(x)
        x=self.activation(x)
        x=self.mc2(x)
        x=self.activation(x)
        x=self.mc3(x)
        return x



def train(model, device, train_loader, loss_criterion, accu_criterion, optimizer, epoch, logwandb):
    batch_time = AverageMeter('Batch', ':.4f')
    data_time = AverageMeter('Data', ':.4f')
    losses = AverageMeter('Loss', ':.4f')
    accues = AverageMeter('Accu', ':.4f')
    progress = ProgressMeter(
        len(train_loader),
        [batch_time, data_time, losses, accues],
        prefix='Epoch: [{}]'.format(epoch))

    model.train()

    end = time.time()

    for i, (data,  target) in enumerate(train_loader):

        data_time.update(time.time() - end)
        data = data.to(device, non_blocking=True)
        target = target.to(device, non_blocking=True)

        output = model(data)

        loss = loss_criterion(output, target)
        losses.update(loss.item(), target.size(0))

        accu = accu_criterion(output, target)
        accues.update(accu.item(), target.size(0))

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        batch_time.update(time.time() - end)
        end = time.time()

        if i % 10 == 0:
            progress.display(i)

        if logwandb:
            wandb.log({'train_loss': losses.avg, 'train_accu': accues.avg})

    return losses.avg


def validate(model, device, test_loader, loss_criterion, accu_criterion, logwandb, test=False, normalizer=None):
    batch_time = AverageMeter('Batch', ':.4f')
    losses = AverageMeter('Loss', ':.4f')
    accues = AverageMeter('Accu', ':.4f')
    progress = ProgressMeter(
        len(test_loader),
        [batch_time, losses, accues],
        prefix='Val: ')

    model.eval()

    with torch.no_grad():
        end = time.time()
        if test:
            y_pred, y_true = [], []
        for i, (data,  target) in enumerate(test_loader):

            data = data.to(device, non_blocking=True)
            target = target.to(device, non_blocking=True)

            output = model(data)

            loss = loss_criterion(output, target)
            losses.update(loss.item(), target.size(0))

            if normalizer is not None:
                output = normalizer.denormalize(output)
                target = normalizer.denormalize(target)

            accu = accu_criterion(output, target)
            accues.update(accu.item(), target.size(0))

            batch_time.update(time.time() - end)
            end = time.time()

            if i % 10 == 0:
                progress.display(i)

            if logwandb:
                wandb.log({'test_loss': losses.avg, 'test_accu': accues.avg})

            if test:
                y_pred.append(output.detach().cpu().numpy())
                y_true.append(target.detach().cpu().numpy())

    if not test:
        return losses.avg, accues.avg
    else:
        return np.vstack(y_pred), np.vstack(y_true)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--batch-size', type=int, default=64, metavar='N',
                        help='input batch size for training')
    parser.add_argument('--train-ratio', type=float, default=0.8,
                        help='trainning ratio')
    parser.add_argument('--start-epoch', default=1, type=int, metavar='N',
                        help='manual epoch number (useful on restarts)')
    parser.add_argument('--epochs', type=int, default=2000, metavar='N',
                        help='number of epochs to train')

    parser.add_argument('--optim', default='SGD', type=str, metavar='Adam',
                        help='choose an optimizer, SGD or Adam, (default:Adam)')
    parser.add_argument('--lr', type=float, default=1e-2, metavar='LR',
                        help='learning rate')
    parser.add_argument('--momentum', default=0.9, type=float, metavar='M',
                        help='momentum')
    parser.add_argument('--weight-decay', default=1e-4, type=float,
                        metavar='W', help='weight decay',
                        dest='weight_decay')
    parser.add_argument('--gamma', type=float, default=0.1, metavar='M',
                        help='Learning rate step gamma')

    parser.add_argument('--seed', type=int, default=4, metavar='S',
                        help='random seed (default: 3)')
    parser.add_argument('--print-freq', default=20, type=int,
                        metavar='N', help='print frequency (default: 20)')

    parser.add_argument('--resume', default=False, action='store_true',
                        help='Resume from checkpoint')
    parser.add_argument('--num-workers', default=8, type=int)
    parser.add_argument('--drop-last', default=False, type=bool)
    parser.add_argument('--pin-memory', default=True, type=bool)
    parser.add_argument('--wandb', default=False, action='store_true')

    args = parser.parse_args()

    best_accu = 1e6
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(args.seed)
    torch.backends.cudnn.benchmark = True

    dataset = NNDataset()
    # normalizer = Normalizer(dataset.target_scaler.mean_, dataset.target_scaler.var_, device)
    
    n_data = len(dataset)
    train_split = int(n_data * args.train_ratio)
    dataset_train, dataset_val = random_split(
        dataset,
        [train_split, len(dataset) - train_split],
        generator=torch.Generator().manual_seed(args.seed)
    )

    train_loader = DataLoader(
        dataset_train,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        drop_last=args.drop_last,
        pin_memory=args.pin_memory,
        generator=torch.Generator().manual_seed(args.seed),

    )

    val_loader = DataLoader(
        dataset_val,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        drop_last=args.drop_last,
        pin_memory=args.pin_memory,

    )

    data0, _ = dataset[0]
    model = NNModel(data0.size(-1))
    model.to(device)

    if args.optim == 'SGD':
        optimizer = optim.SGD(model.parameters(), lr=args.lr, momentum=args.momentum, weight_decay=args.weight_decay)
    elif args.optim == 'Adam':
        optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    else:
        raise ValueError('Optimizer must be SGD or Adam.')
    loss_criterion = torch.nn.MSELoss()
    accu_criterion = torch.nn.L1Loss()

    if args.resume:
        print("=> loading checkpoint")
        checkpoint = torch.load('checkpoint.pth.tar')
        args.start_epoch = checkpoint['epoch'] + 1
        best_accu = checkpoint['best_accu']
        model.load_state_dict(checkpoint['state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        print("=> loaded checkpoint (epoch {})".format(checkpoint['epoch']))

    if args.wandb:
        wandb.init(project="mlff-project", entity="wayne833")
        wandb.watch(model, log='all', log_freq=100)

    for epoch in range(args.start_epoch, args.epochs + 1):
        train(model, device, train_loader, loss_criterion, accu_criterion, optimizer, epoch, args.wandb)
        val_loss, val_accu = validate(model, device, val_loader, loss_criterion, accu_criterion, args.wandb)
        # scheduler.step()

        is_best = val_accu < best_accu
        if is_best:
            best_accu = val_accu

        if args.wandb:
            wandb.log({'best_accu': best_accu})

        save_checkpoint({
            'epoch': epoch,
            'state_dict': model.state_dict(),
            'best_accu': best_accu,
            'optimizer': optimizer.state_dict(),
            'args': vars(args),
        }, is_best)

    ckpt = torch.load('model_best.pth.tar')
    model.load_state_dict(ckpt['state_dict'])
    model.eval()
    pred, target = validate(model, device, val_loader, loss_criterion, accu_criterion, logwandb=False, test=True)
    with open('pred_target.dat', 'w') as f:
        for row in range(len(pred)):
            for j in range(len(pred[row])):
                f.write("%9.6f  %9.6f\n" % (pred[row][j],target[row][j]));



def save_checkpoint(state, is_best, filename='checkpoint.pth.tar'):
    torch.save(state, filename)
    if is_best:
        shutil.copyfile(filename, 'model_best.pth.tar')


if __name__ == '__main__':
    main()
