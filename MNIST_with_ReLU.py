# %% imports
# pytorch
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from torchvision import transforms, datasets
from torch.utils.data import Dataset, DataLoader
import numpy as np

# pyplot
import matplotlib.pyplot as plt


# %% Noisy MNIST dataset
class Noisy_MNIST(Dataset):
    # initialization of the dataset
    def __init__(self, split, data_loc, noise=0.5):
        # save the input parameters
        self.split = split
        self.data_loc = data_loc
        self.noise = noise

        if self.split == 'train':
            Train = True
        else:
            Train = False

        # get the original MNIST dataset
        Clean_MNIST = datasets.MNIST(self.data_loc, train=Train, download=True)

        # reshuffle the test set to have digits 0-9 at the start
        if self.split == 'train':
            data = Clean_MNIST.data.unsqueeze(1)
        else:
            data = Clean_MNIST.data.unsqueeze(1)
            idx = torch.load('test_idx.tar')
            data[:, :] = data[idx, :]

        # reshape and normalize
        resizer = transforms.Resize(32)
        resized_data = resizer(data) * 1.0
        normalized_data = 2 * (resized_data / 255) - 1
        # normalized_data = (resized_data - 33)/74

        # create the data
        self.Clean_Images = normalized_data
        self.Noisy_Images = normalized_data + torch.randn(normalized_data.size()) * self.noise
        self.Labels = Clean_MNIST.targets

    # return the number of examples in this dataset
    def __len__(self):
        return self.Labels.size(0)

    # create a a method that retrieves a single item form the dataset
    def __getitem__(self, idx):
        clean_image = self.Clean_Images[idx, :, :, :]
        noisy_image = self.Noisy_Images[idx, :, :, :]
        label = self.Labels[idx]

        return clean_image, noisy_image, label


class FCN(nn.Module):
    def __init__(self):
        super(FCN, self).__init__()
        self.out_1 = nn.Linear(1024, 2000)
        self.out_2 = nn.Linear(2000, 512)
        self.out_3 = nn.Linear(512, 256)
        self.out_4 = nn.Linear(256, 1024)
        # self.out_5 = nn.Linear(64, 10)

    ##use drop out and batch norm
    # self.dropout=nn.Dropout(0.5)
    #     self.bn1=nn.BatchNorm1d(2000)
    #     self.bn2 = nn.BatchNorm1d(512)
    #     self.bn3 = nn.BatchNorm1d(256)
        #self.bn4 = nn.BatchNorm1d(1024)
    # self.softmax=nn.LogSoftmax(dim=1)

    ##activation function
    def Relu(self,x):
        m=x.data
        m[m[:,:]<0]= 0
        x.data=m
        return x

    def forward(self, x):
        x = x.view(-1, 1024)
        # output1 = F.relu(self.bn1(self.out_1(x)))
        # output1=self.dropout(output1)
        # output2=  F.relu(self.bn2(self.out_2(output1)))
        # output2 = self.dropout(output2)
        # output3=  F.relu(self.bn3(self.out_3(output2)))
        # output3=self.dropout(output3)
        # output4 = F.relu(self.bn4(self.out_4(output3)))
        # output4 = self.dropout(output4)
        # output5=self.out_5(output4)
        # output5=self.softmax(output5)
        output1 = self.Relu((self.out_1(x)))
        #output1[:,:]=0
        #print(output1)
        output2 = (self.out_2(output1))
        output3 = (self.out_3(output2))
        output4 = self.out_4(output3)
        #print(output4)
        # output5 = self.out_5(output4)
        # output5 = self.softmax(output4)
        return output4


fcn = FCN()
# optimizer = torch.optim.SGD(fcn.parameters(),
#                             lr=0.002,
#                             momentum=0.9)
optimizer = torch.optim.Adam(fcn.parameters(), lr=0.0001)

# lossfn = F.mse_loss()


# %% dataloader for the Noisy MNIST dataset
def create_dataloaders(data_loc, batch_size):
    Noisy_MNIST_train = Noisy_MNIST("train", data_loc)
    Noisy_MNIST_test = Noisy_MNIST("test", data_loc)

    train_size = int(0.8 * len(Noisy_MNIST_train))
    validation_size = len(Noisy_MNIST_train) - train_size
    Noisy_MNIST_train, Noisy_MNIST_validation = torch.utils.data.random_split(Noisy_MNIST_train,
                                                                              [train_size, validation_size])

    Noisy_MNIST_train_loader = DataLoader(Noisy_MNIST_train, batch_size=batch_size, shuffle=True, drop_last=False)
    Noisy_MNIST_test_loader = DataLoader(Noisy_MNIST_test, batch_size=batch_size, shuffle=False, drop_last=False)
    Noisy_MNIST_validation_loader = DataLoader(Noisy_MNIST_validation, batch_size=batch_size, shuffle=True,
                                               drop_last=False)

    return Noisy_MNIST_train_loader, Noisy_MNIST_test_loader, Noisy_MNIST_validation_loader


# %% test if the dataloaders work
if __name__ == "__main__":
    # define parameters
    data_loc = 'D://5LSL0/dataset'  # change the datalocation to something that works for you
    batch_size = 64

    # get dataloader
    train_loader, test_loader, validation_loader = create_dataloaders(data_loc, batch_size)

    # train the network
    EPOCH = 20
    t = 0
    trainloss = []
    valiloss = []
    for epoch in range(EPOCH):
        ##save untrained model
        if epoch == 0:
            torch.save(fcn, './trained_model/fcnuntrained.pth')

        loss = 0
        for t, (xclean, xnoise, y) in enumerate(train_loader):
            xn = xnoise.to(device=torch.device('cpu'), dtype=torch.float32)
            xc = xclean.to(device=torch.device('cpu'), dtype=torch.float32)

            xn = Variable(xn)
            xc = Variable(xc)

            ##train the data
            outputtrain = fcn(xn)
            #print((outputtrain))
            loss = F.mse_loss(outputtrain, xc.view(-1, 32 * 32))  # calculate loss
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # print training loss
        print('now epoch :  ', epoch, '   |  trainloss : %.4f ' % loss.item())
        ##save trainloss to list
        trainloss.append(loss.item())

        ##validation loss
        for v, (vxclean, vxnoise, vy) in enumerate(validation_loader):
            vxn = vxnoise.to(device=torch.device('cpu'), dtype=torch.float32)  # move to device, e.g. GPU
            vxc = vxclean.to(device=torch.device('cpu'), dtype=torch.float32)
            vxn = Variable(vxn)
            vxc = Variable(vxc)
            ##test loss in validation set
            outvalidation = fcn(vxn)
            validationloss = F.mse_loss(outvalidation, vxc.view(-1, 32 * 32))
        print('now epoch :  ', epoch, '   |  validationloss : %.4f ' % validationloss.item())
        ##save the validation loss to list
        valiloss.append(validationloss.item())
    # print(len(trainloss))
    ##save model after training
    torch.save(fcn, './trained_model/fcntrained.pth')

    ##load the model for test
    modeltestuntrain = torch.load('./trained_model/fcnuntrained.pth')
    modeltesttrain = torch.load('./trained_model/fcntrained.pth')
    testlossuntrain = []
    testlosstrain = []

    # calculate test loss
    for test, (txclean, txnoise, ty) in enumerate(test_loader):
        txn = txnoise.to(device=torch.device('cpu'), dtype=torch.float32)  # move to device, e.g. GPU
        txc = txclean.to(device=torch.device('cpu'), dtype=torch.float32)
        txn = Variable(txn)
        txc = Variable(txc)
        ##test loss in untrain model
        testoutuntrain = modeltestuntrain(txn)
        testlossuntrain.append(F.mse_loss(testoutuntrain, txc.view(-1, 32 * 32)).item())

        ##test loss in trained model
        testouttrain = modeltesttrain(txn)
        testlosstrain.append(F.mse_loss(testouttrain, txc.view(-1, 32 * 32)).item())
        # y_pred = torch.max(testout, 1)[1].data.squeeze()
        # accuracy = sum(txc == y_pred).item() / txc.size(0)

    # print('testlossbefore training : %.4f ' % testlossuntrain.item())
    # print('testloss after training: %.4f ' % testlosstrain.item())

    ##plot training loss and validation loss
    fig = plt.figure(figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.subplot(2, 1, 1)
    plt.plot(range(EPOCH), trainloss,
             label='training')
    plt.plot(range(EPOCH), valiloss,
             label='validation', linestyle='--')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(range(len(testlossuntrain)), testlossuntrain,
             label='testsetuntrain')
    plt.plot(range(len(testlosstrain)), testlosstrain,
             label='testsettrained', linestyle='--')
    plt.xlabel('testdata')
    plt.ylabel('Loss')
    plt.legend()

    plt.show()

    # get some examples
    examples = enumerate(test_loader)
    _, (x_clean_example, x_noisy_example, labels_example) = next(examples)
    # use these example images througout the assignment as the first 10 correspond to the digits 0-9

    # show the examples in a plot
    plt.figure(figsize=(12, 3))
    for i in range(10):
        if i == 4:
            plt.subplot(4, 10, i + 1, title='clean test data')
        else:
            plt.subplot(4, 10, i + 1)
        plt.imshow(x_clean_example[i, 0, :, :], cmap='gray')
        plt.xticks([])
        plt.yticks([])

        if i == 4:
            plt.subplot(4, 10, i + 11, title='noise test data')
        else:
            plt.subplot(4, 10, i + 11)
        plt.imshow(x_noisy_example[i, 0, :, :], cmap='gray')
        plt.xticks([])
        plt.yticks([])

        ##print test set result from untrained model
        outuntrain = modeltestuntrain(x_noisy_example[i, 0, :, :])
        outuntrain = outuntrain.reshape(32, 32)
        # print(x_noisy_example[i,0,:,:],outuntrain.data)
        if i == 4:
            plt.subplot(4, 10, i + 21, title='  test result from untrained model')
        else:
            plt.subplot(4, 10, i + 21)
        plt.imshow(outuntrain.data, cmap='gray')
        plt.xticks([])
        plt.yticks([])

        ##print test set result from trained model
        outtrain = modeltesttrain(x_noisy_example[i, 0, :, :])
        outtrain = outtrain.reshape(32, 32)
        # print(x_noisy_example[i,0,:,:],outuntrain.data)
        if i == 4:
            plt.subplot(4, 10, i + 31, title='  test result from trained model')
        else:
            plt.subplot(4, 10, i + 31)
        plt.imshow(outtrain.data, cmap='gray')
        plt.xticks([])
        plt.yticks([])

    plt.tight_layout()
    plt.savefig("data_examples.png", dpi=300, bbox_inches='tight')
    plt.show()
