few_shot_examples = """
    Task: Implement a class to process and split genomic data.
    Input: path_to_data
    Code:
    ```python
    import pandas as pd
    import numpy as np

    class DataProcessor:
        def __init__(self, path_to_data):
            self.path = path_to_data

        def concat_data(self):
            cell_types_v = ['GM', 'H1', 'K562', 'MCF7']
            positive, type_1_negative, type_2_negative, type_3_negative = [], [], [], []

            for cell_type in cell_types_v:
                positive.append(pd.read_csv(self.path + cell_type + '_insulator_pos_withCTCF.fa', sep=">chr*", header=None, engine='python').values[1::2][:,0])
                type_1_negative.append(pd.read_csv(self.path + cell_type + '_type1.fa', sep=">chr*", header=None, engine='python').values[1::2][:,0])
                type_2_negative.append(pd.read_csv(self.path + cell_type + '_type2.fa', sep=">chr*", header=None, engine='python').values[1::2][:,0])
                type_3_negative.append(pd.read_csv(self.path + cell_type + '_type3.fa', sep=">chr*", header=None, engine='python').values[1::2][:,0])

            return positive, type_1_negative, type_2_negative, type_3_negative

        def split(self, file, size=0.1):
            len_v = int(len(file) * size)
            np.random.seed(42)
            np.random.shuffle(file)
            train, test = file[len_v:], file[:len_v]
            train, val = train[len_v:], train[:len_v]
            return train, test, val
    ```

    Task: Implement a function to compute the reverse complement of a DNA sequence.
    Input: 'ATGC'
    Output: 'GCAT'
    Code:
    ```python
    def RC(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    t = ''
    for base in seq:
        t = complement[base] + t
    return t
    ```

    Task: Implement a PyTorch Dataset class for siamese neural network data.
    Input: data, label
    Code:
    ```python
    import torch
    from torch.utils.data import Dataset

    class Data_siam(Dataset):
        def __init__(self, data, label):
            self.data = data
            self.label = label

        def __len__(self):
            return len(self.data)

        def __getitem__(self, index):
            seq = self.data[index]
            rc_seq = RC(seq)
            ctr = 0
            ar1 = np.zeros((2000, 4))
            for base in seq:
                if base == 'A' or base == 'a':
                    ar1[ctr, 0] = 1
                elif base == 'T' or base == 't':
                    ar1[ctr, 1] = 1
                elif base == 'C' or base == 'c':
                    ar1[ctr, 2] = 1
                elif base == 'G' or base == 'g':
                    ar1[ctr, 3] = 1
                ctr += 1

            ar2 = np.zeros((2000, 4))
            ctr = 0
            for base in rc_seq:
                if base == 'A' or base == 'a':
                    ar2[ctr, 0] = 1
                elif base == 'T' or base == 't':
                    ar2[ctr, 1] = 1
                elif base == 'C' or base == 'c':
                    ar2[ctr, 2] = 1
                elif base == 'G' or base == 'g':
                    ar2[ctr, 3] = 1
                ctr += 1

            ar1 = torch.tensor(ar1).float().permute(1, 0)
            ar2 = torch.tensor(ar2).float().permute(1, 0)
            label = torch.tensor(self.label).float()

            return ar1, ar2, label
    ```

    Task: Implement a neural network class with LSTM in PyTorch.
    Input: N/A
    Code:
    ```python
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    class NNetwork_wLSTM(nn.Module):
        def __init__(self):
            super(NNetwork_wLSTM, self).__init__()
            self.Conv1 = nn.Conv1d(in_channels=4, out_channels=160, kernel_size=31)
            self.Maxpool1 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv2 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=20)
            self.Maxpool2 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv3 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=6)
            self.Maxpool3 = nn.MaxPool1d(kernel_size=8, stride=6)
            self.BiLSTM = nn.LSTM(input_size=160, hidden_size=160, num_layers=2,
                                batch_first=True, dropout=0.5, bidirectional=True)
            self.Drop1 = nn.Dropout(p=0.3)
            self.Linear1 = nn.Linear(79*320, 925)
            self.Linear2 = nn.Linear(925, 925)
            self.Linear3 = nn.Linear(925, 1)

        def forward_one(self, input):
            x = self.Conv1(input)
            x = F.relu(x)
            x = self.Maxpool1(x)
            x = self.Conv2(x)
            x = F.relu(x)
            x = self.Maxpool2(x)
            x = self.Conv3(x)
            x = F.relu(x)
            x = self.Maxpool3(x)
            x_x = torch.transpose(x, 1, 2)
            x, (h_n, h_c) = self.BiLSTM(x_x)
            x = x.contiguous().view(-1, 79*320)
            x = self.Drop1(x)
            x = self.Linear1(x)
            x = F.relu(x)
            x = self.Drop1(x)
            x = self.Linear2(x)
            x = F.relu(x)
            x = self.Linear3(x)
            return x

        def forward(self, x1, x2):
            out1 = self.forward_one(x1)
            out2 = self.forward_one(x2)
            out = (out1 + out2) / 2
            return torch.sigmoid(out)
    ```

    Task: Implement a class for a simple convolutional neural network in PyTorch.
    Input: N/A
    Code:
    ```python
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    class NNetwork(nn.Module):
        def __init__(self):
            super(NNetwork, self).__init__()
            self.Conv1 = nn.Conv1d(in_channels=4, out_channels=160, kernel_size=31)
            self.Maxpool1 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv2 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=20)
            self.Maxpool2 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv3 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=6)
            self.Maxpool3 = nn.MaxPool1d(kernel_size=8, stride=6)
            self.Drop1 = nn.Dropout(p=0.3)
            self.Linear1 = nn.Linear(79*160, 925)
            self.Linear2 = nn.Linear(925, 925)
            self.Linear3 = nn.Linear(925, 1)

        def forward_one(self, input):
            x = self.Conv1(input)
            x = F.relu(x)
            x = self.Maxpool1(x)
            x = self.Conv2(x)
            x = F.relu(x)
            x = self.Maxpool2(x)
            x = self.Conv3(x)
            x = F.relu(x)
            x = self.Maxpool3(x)
            x = torch.flatten(x, 1)
            x = self.Drop1(x)
            x = self.Linear1(x)
            x = F.relu(x)
            x = self.Drop1(x)
            x = self.Linear2(x)
            x = F.relu(x)
            x = self.Linear3(x)
            return x

        def forward(self, x1, x2):
            out1 = self.forward_one(x1)
            out2 = self.forward_one(x2)
            out = (out1 + out2) / 2
            return torch.sigmoid(out)
    ```

    Task: Implement a function to calculate precision, recall, and F1-score from predictions and ground truth labels.
    Input: y_true, y_pred
    Code:
    ```python
    from sklearn.metrics import precision_score, recall_score, f1_score

    def compute_metrics(y_true, y_pred):
        precision = precision_score(y_true, y_pred)
        recall = recall_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)
        return precision, recall, f1
    ```
     Task: Implement a convolutional neural network in PyTorch.
    Input: None
    Code:
    ```python
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    class NNetwork(nn.Module):
        def __init__(self):
            super(NNetwork, self).__init__()
            self.Conv1 = nn.Conv1d(in_channels=4, out_channels=160, kernel_size=31)
            self.Maxpool1 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv2 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=20)
            self.Maxpool2 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv3 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=6)
            self.Maxpool3 = nn.MaxPool1d(kernel_size=8, stride=6)
            self.Drop1 = nn.Dropout(p=0.3)
            self.Linear1 = nn.Linear(79*160, 925)
            self.Linear2 = nn.Linear(925, 925)
            self.Linear3 = nn.Linear(925, 1)

        def forward_one(self, input):
            x = self.Conv1(input)
            x = F.relu(x)
            x = self.Maxpool1(x)
            x = self.Conv2(x)
            x = F.relu(x)
            x = self.Maxpool2(x)
            x = self.Conv3(x)
            x = F.relu(x)
            x = self.Maxpool3(x)
            x = torch.flatten(x, 1)
            x = self.Drop1(x)
            x = self.Linear1(x)
            x = F.relu(x)
            x = self.Drop1(x)
            x = self.Linear2(x)
            x = F.relu(x)
            x = self.Linear3(x)
            return x

        def forward(self, x1, x2):
            out1 = self.forward_one(x1)
            out2 = self.forward_one(x2)
            out = (out1 + out2) / 2
            return torch.sigmoid(out)
    ```

    Task: Count the number of parameters in a PyTorch model.
    Input: NNetwork()
    Output: Number of Parameters - NNetwork :  170674
    Code:
    ```python
    def count_parameters(model):
        return sum(p.numel() for p in model.parameters() if p.requires_grad)

    print('Number of Parameters - NNetwork : ', count_parameters(NNetwork()))
    ```

    Task: Compute the precision-recall AUC score from true labels and predicted probabilities.
    Input: y_true, y_proba
    Output: Precision-Recall AUC score
    Code:
    ```python
    def prcs(y, y_proba):
        from sklearn.metrics import precision_recall_curve
        from sklearn.metrics import auc
        lr_precision, lr_recall, _ = precision_recall_curve(y, y_proba)
        lr_auc = auc(lr_recall, lr_precision)
        return lr_auc
    ```

    Task: Compute the ROC AUC score from true labels and predicted probabilities.
    Input: y_true, y_proba
    Output: ROC AUC score
    Code:
    ```python
    def rocs(y, y_proba):
        from sklearn.metrics import roc_auc_score
        lr_auc = roc_auc_score(y, y_proba)
        return lr_auc
    ```

    Task: Implement a PyTorch Trainer class for model training and evaluation.
    Input: train_data, val_data, model, num_epochs, batch_size, learning_rate, weight_decay, pretrain_path, model_path
    Code:
    ```python
    import torch
    import torch.nn as nn
    from tqdm import tqdm
    from sklearn.metrics import classification_report

    class Trainer:
        def __init__(self, train_data, val_data, model, num_epochs, batch_size,
                     learning_rate, weight_decay, pretrain_path, model_path):
            self.train_data = train_data
            self.val_data = val_data
            self.model_path = model_path
            self.num_epochs = num_epochs
            self.batch_size = batch_size
            self.learning_rate = learning_rate
            self.weight_decay = weight_decay
            self.pretrain_path = pretrain_path
            self.model = model

        def train(self):
            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            train_loader = torch.utils.data.DataLoader(self.train_data, shuffle=True, batch_size=self.batch_size)
            val_loader = torch.utils.data.DataLoader(self.val_data, shuffle=True, batch_size=self.batch_size)
            model = self.model.to(device)
            min_loss = 100
            optimizer = torch.optim.Adam(model.parameters(), lr=self.learning_rate)
            criterion = nn.BCELoss()

            for epoch in range(1, self.num_epochs + 1):
                print("Epoch {}".format(epoch))
                model.train()
                train_acc = 0
                train_loss = 0
                for data, label in tqdm(train_loader):
                    data, label = data.to(device), label.to(device)
                    output = model.forward_one(data).squeeze()
                    loss = criterion(torch.sigmoid(output), label)
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                    train_loss += loss.item()
                    y_pred = (output > 0.5).float()
                    train_acc += torch.sum(y_pred == label)

                loss = train_loss / len(train_loader)
                accuracy = int(train_acc / (len(train_loader.dataset)) * 100)
                print('\n Train Data: Average Train Loss: {:.4f}, Train Accuracy: {}/{} ({}%)'.format(loss, train_acc, len(train_loader.dataset), accuracy))

                y_true = []
                y_proba = []
                y_pred = []
                model.eval()
                val_loss = 0
                val_accuracy = 0

                with torch.no_grad():
                    for data, target in tqdm(val_loader):
                        data, target = data.to(device), target.to(device)
                        output = model.forward_one(data).squeeze()
                        y_hat = (output).cpu().numpy()
                        loss = criterion(torch.sigmoid(output), target)
                        val_loss += loss.item() * data.size(0)
                        y_pred_ = (output > 0.5).float()
                        val_accuracy += sum(y_pred_ == target)
                        for i in range(len(y_pred_)):
                            y_true.append(float(target[i]))
                            y_pred.append(float(y_pred_[i]))
                            y_proba.append(float(y_hat[i]))

                    loss = val_loss / len(val_loader.dataset)
                    accuracy = val_accuracy / len(val_loader.dataset)
                    prc = prcs(y_true, y_proba)
                    roc = rocs(y_true, y_proba)
                    print('Validation -> AUPRC : {:.4f} , AUROC : {:.4f} , Loss : {:.4f}'.format(prc, roc, loss))
                    print('#')
                    print(classification_report(y_true, y_pred, target_names=['Negative', 'Positive']))

                if loss < min_loss:
                    min_loss = loss
                    torch.save(model.state_dict(), self.model_path)
    ```
    Task: Implement a `NNetwork_wLSTM` class with convolutional and LSTM layers in PyTorch.
    Input: None
    Code:
    ```python
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    class NNetwork_wLSTM(nn.Module):
        def __init__(self):
            super(NNetwork_wLSTM, self).__init__()
            self.Conv1 = nn.Conv1d(in_channels=4, out_channels=160, kernel_size=31)
            self.Maxpool1 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv2 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=20)
            self.Maxpool2 = nn.MaxPool1d(kernel_size=2, stride=2)
            self.Conv3 = nn.Conv1d(in_channels=160, out_channels=160, kernel_size=6)
            self.Maxpool3 = nn.MaxPool1d(kernel_size=8, stride=6)
            self.Drop1 = nn.Dropout(p=0.3)

            self.LSTM = nn.LSTM(input_size=160, hidden_size=128, num_layers=2, bidirectional=True, batch_first=True)

            self.Linear1 = nn.Linear(256, 925)
            self.Linear2 = nn.Linear(925, 925)
            self.Linear3 = nn.Linear(925, 1)

        def forward(self, x):
            x = self.Conv1(x)
            x = F.relu(x)
            x = self.Maxpool1(x)
            x = self.Conv2(x)
            x = F.relu(x)
            x = self.Maxpool2(x)
            x = self.Conv3(x)
            x = F.relu(x)
            x = self.Maxpool3(x)
            x = torch.flatten(x, 1)
            x = self.Drop1(x)
            x, _ = self.LSTM(x.unsqueeze(1))
            x = x[:, -1, :]
            x = self.Linear1(x)
            x = F.relu(x)
            x = self.Drop1(x)
            x = self.Linear2(x)
            x = F.relu(x)
            x = self.Linear3(x)
            return torch.sigmoid(x)
    ```

    Task: Implement a class to calculate the classification report for model evaluation.
    Input: y_true, y_pred
    Output: Classification report
    Code:
    ```python
    from sklearn.metrics import classification_report

    def get_classification_report(y_true, y_pred):
        return classification_report(y_true, y_pred, target_names=['Negative', 'Positive'])
    ```

    Task: Compute the average loss and accuracy of a PyTorch model over a dataset.
    Input: model, data_loader, criterion
    Output: Average loss and accuracy
    Code:
    ```python
    def evaluate_model(model, data_loader, criterion, device):
        model.eval()
        total_loss = 0
        correct = 0
        with torch.no_grad():
            for data, target in data_loader:
                data, target = data.to(device), target.to(device)
                output = model(data).squeeze()
                loss = criterion(torch.sigmoid(output), target)
                total_loss += loss.item() * data.size(0)
                pred = (output > 0.5).float()
                correct += pred.eq(target).sum().item()
        
        avg_loss = total_loss / len(data_loader.dataset)
        accuracy = correct / len(data_loader.dataset)
        return avg_loss, accuracy
    ```

    Task: Save and load a PyTorch model's state dictionary.
    Input: model, model_path
    Code:
    ```python
    def save_model(model, model_path):
        torch.save(model.state_dict(), model_path)

    def load_model(model, model_path):
        model.load_state_dict(torch.load(model_path))
        model.eval()
        return model
    ```

    Task: Train a PyTorch model with early stopping based on validation loss.
    Input: Trainer class instance
    Code:
    ```python
    class TrainerWithEarlyStopping(Trainer):
        def __init__(self, train_data, val_data, model, num_epochs, batch_size,
                     learning_rate, weight_decay, pretrain_path, model_path, patience):
            super().__init__(train_data, val_data, model, num_epochs, batch_size, learning_rate, weight_decay, pretrain_path, model_path)
            self.patience = patience
            self.best_loss = float('inf')
            self.patience_counter = 0

        def train(self):
            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            train_loader = torch.utils.data.DataLoader(self.train_data, shuffle=True, batch_size=self.batch_size)
            val_loader = torch.utils.data.DataLoader(self.val_data, shuffle=True, batch_size=self.batch_size)
            model = self.model.to(device)
            optimizer = torch.optim.Adam(model.parameters(), lr=self.learning_rate)
            criterion = nn.BCELoss()

            for epoch in range(1, self.num_epochs + 1):
                print("Epoch {}".format(epoch))
                model.train()
                train_acc = 0
                train_loss = 0
                for data, label in tqdm(train_loader):
                    data, label = data.to(device), label.to(device)
                    output = model.forward_one(data).squeeze()
                    loss = criterion(torch.sigmoid(output), label)
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                    train_loss += loss.item()
                    y_pred = (output > 0.5).float()
                    train_acc += torch.sum(y_pred == label)

                loss = train_loss / len(train_loader)
                accuracy = int(train_acc / (len(train_loader.dataset)) * 100)
                print('\n Train Data: Average Train Loss: {:.4f}, Train Accuracy: {}/{} ({}%)'.format(loss, train_acc, len(train_loader.dataset), accuracy))

                val_loss, val_accuracy = evaluate_model(model, val_loader, criterion, device)
                print('Validation -> Loss : {:.4f}, Accuracy: {:.4f}'.format(val_loss, val_accuracy))

                if val_loss < self.best_loss:
                    self.best_loss = val_loss
                    self.patience_counter = 0
                    save_model(model, self.model_path)
                else:
                    self.patience_counter += 1
                    if self.patience_counter >= self.patience:
                        print("Early stopping")
                        break
    ```

        Task: Flatten nested lists and split data for training, testing, and validation.
    Input: 
    ```python
    def flatten(t):
        return [item for sublist in t for item in sublist]

    ptrain,ptest,pval = process.split(flatten(positive))
    type_1train,type_1test,type_1val = process.split(flatten(type_1))
    type_2train,type_2test,type_2val = process.split(flatten(type_2))
    type_3train,type_3test,type_3val = process.split(flatten(type_3))
    ```
    Output: 
    ```python
    # Nested lists are flattened and split into ptrain, ptest, pval, type_1, type_2, and type_3 datasets.
    ```

    Task: Create PyTorch DataLoader for training and validation datasets.
    Input: 
    ```python
    Ptrain_data = Data(ptrain,1)
    Pval_data = Data(pval,1)

    type1_train_data = Data(type_1train,0)
    type1_val_data = Data(type_1val,0)

    type2_train_data = Data(type_2train,0)
    type2_val_data = Data(type_2val,0)

    type3_train_data = Data(type_3train,0)
    type3_val_data = Data(type_3val,0)

    train_data_loader = torch.utils.data.ConcatDataset([Ptrain_data,type1_train_data,type2_train_data,type3_train_data])
    val_data_loader = torch.utils.data.ConcatDataset([Pval_data,type1_val_data, type2_val_data, type3_val_data])
    ```
    Output: 
    ```python
    # DataLoader objects are created for training and validation datasets.
    ```

    Task: Train a PyTorch model using the Trainer class.
    Input: 
    ```python
    pre_train_path = ''
    model_path = '/content/drive/MyDrive/model.pt'
    model = NNetwork()

    epochs = 15
    batch_size = 64
    lr = 5e-5
    weight_decay = 5e-4

    trainer = Trainer(train_data_loader, val_data_loader, model,
                      epochs, batch_size, lr, weight_decay, pre_train_path, model_path)

    trainer.train()
    ```
    Output: 
    ```python
    # The model is trained with specified parameters and the best model is saved to the given path.
    ```

    Task: Test the PyTorch model and evaluate its performance.
    Input: 
    ```python
    class Tester():
        def __init__(self, model, model_weight_path, test_data, batch_size):
            self.model = model
            self.model_weight_path = model_weight_path
            self.test_data = test_data
            self.batch_size = batch_size

        def test(self):
            test_loader = torch.utils.data.DataLoader(self.test_data, batch_size=self.batch_size, shuffle=True)
            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            model = self.model.to(device)
            model.load_state_dict(torch.load(self.model_weight_path))

            test_loss = 0
            correct = 0
            y_pred = np.array([])
            y_true = np.array([])
            y_proba = np.array([])
            model.eval()

            with torch.no_grad():
                for X1, X2, y in tqdm(test_loader):
                    X1, X2, y = X1.to(device), X2.to(device), y.to(device)
                    output = model(X1.float(), X2.float())
                    y_hat = output
                    y = y.float()

                    y_true = np.concatenate((y_true, y.cpu().numpy()))
                    y_pred = np.concatenate((y_pred, y_hat.cpu().numpy().flatten()))
            return y_true, y_pred
    ```
    Output: 
    ```python
    # The model is tested on the given test data and predictions are returned.
    ```

    Task: Plot Precision-Recall and ROC curves for model evaluation.
    Input: 
    ```python
    class plot_prc_roc():
        def __init__(self, y_true, y_pred, y_true1, y_pred1, y_true2, y_pred2, y_true3, y_pred3, out_path):
            self.y_true = y_true
            self.y_true1 = y_true1
            self.y_true2 = y_true2
            self.y_true3 = y_true3
            self.y_pred = y_pred
            self.y_pred1 = y_pred1
            self.y_pred2 = y_pred2
            self.y_pred3 = y_pred3
            self.out_path = out_path

        def plot_prc(self):
            pyplot.figure(dpi=500)

            lr_precision, lr_recall, _ = precision_recall_curve(self.y_true, self.y_pred)
            lr_auco = auc(lr_recall, lr_precision)
            pyplot.plot(lr_recall, lr_precision)

            lr_precision, lr_recall, _ = precision_recall_curve(self.y_true1, self.y_pred1)
            lr_auc1 = auc(lr_recall, lr_precision)
            pyplot.plot(lr_recall, lr_precision)

            lr_precision, lr_recall, _ = precision_recall_curve(self.y_true2, self.y_pred2)
            lr_auc2 = auc(lr_recall, lr_precision)
            pyplot.plot(lr_recall, lr_precision)

            lr_precision, lr_recall, _ = precision_recall_curve(self.y_true3, self.y_pred3)
            lr_auc3 = auc(lr_recall, lr_precision)
            pyplot.plot(lr_recall, lr_precision)

            pyplot.xlabel('Recall', fontsize = 13)
            pyplot.ylabel('Precision', fontsize = 13)
            pyplot.title('Precision-Recall curve',fontsize = 15)
            pyplot.legend(['Overall (AUC = %0.3f)' % lr_auco, 'Non-anchor type1 (AUC = %0.3f)' % lr_auc1,'Non-anchor type2 (AUC = %0.3f)' % lr_auc2, 'Non-anchor type3 (AUC = %0.3f)' % lr_auc3 ])
            pyplot.show()

        def plot_roc(self):
            pyplot.figure(dpi=500)

            lr_auco = roc_auc_score(self.y_true, self.y_pred)
            lr_fpr, lr_tpr, _ = roc_curve(self.y_true, self.y_pred)
            pyplot.plot(lr_fpr, lr_tpr)

            lr_auc1 = roc_auc_score(self.y_true1, self.y_pred1)
            lr_fpr, lr_tpr, _ = roc_curve(self.y_true1, self.y_pred1)
            pyplot.plot(lr_fpr, lr_tpr)

            lr_auc2 = roc_auc_score(self.y_true2, self.y_pred2)
            lr_fpr, lr_tpr, _ = roc_curve(self.y_true2, self.y_pred2)
            pyplot.plot(lr_fpr, lr_tpr)

            lr_auc3 = roc_auc_score(self.y_true3, self.y_pred3)
            lr_fpr, lr_tpr, _ = roc_curve(self.y_true3, self.y_pred3)
            pyplot.plot(lr_fpr, lr_tpr)

            pyplot.plot([0, 1], [0, 1], color = 'black', linewidth = 1, linestyle = 'dashed')

            pyplot.xlabel('False Positive Rate',fontsize = 13)
            pyplot.ylabel('True Positive Rate',fontsize = 13)
            pyplot.title('Receiver Operating Characteristic curve', fontsize = 15)
            pyplot.legend(['Overall (AUC = %0.3f)' % lr_auco, 'Non-anchor type1 (AUC = %0.3f)' % lr_auc1,'Non-anchor type2 (AUC = %0.3f)' % lr_auc2, 'Non-anchor type3 (AUC = %0.3f)' % lr_auc3 ])
            pyplot.show()
    ```
    Output: 
    ```python
    # Precision-Recall and ROC curves are plotted to evaluate model performance.
    ```
"""