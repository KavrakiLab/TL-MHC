from torch import nn
import torch.nn.functional as F

class Pretrained_BAEL(nn.Module):
    def __init__(self, num_features, latent_dim):
        super().__init__()

        self.fc1 = nn.Linear(num_features, latent_dim * 2)
        self.fc2 = nn.Linear(latent_dim * 2, latent_dim)
        self.fc3 = nn.Linear(latent_dim, 2)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x