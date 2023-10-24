import numpy as np


class MyGaussianNBClassifier:
    def __init__(self, priors=None):
        self.classes = None
        self.priors = priors
        self.mu = None
        self.sigma = None
        pass

    def fit(self, X, y):
        self.classes, n_of_occurrences = np.unique(y, return_counts=True)
        self.mu = np.zeros((self.classes.shape[0], X.shape[1]))
        self.sigma = np.zeros((self.classes.shape[0], X.shape[1]))
        if self.priors is None:
            self.priors = np.array([float (n_of_occur) / y.shape[0] for n_of_occur in n_of_occurrences])

        x_for_each_class = [[] for _ in self.classes]

        for x_, y_ in zip(X,y):
            x_for_each_class[self.classes[y_]].append(x_)

        for i, xs in enumerate(x_for_each_class):
            # test = np.mean(np.array(xs), axis=0)
            self.mu[i, :] = np.mean(np.array(xs), axis=0)
            self.sigma[i, :] = np.sqrt(np.var(np.array(xs), axis=0))

    def predict(self, X):
        return np.array([self.classes[ind] for ind in np.argmax(self.predict_proba(X), axis=1)])

    def predict_proba(self, X):
        proba = np.zeros((X.shape[0], self.classes.shape[0]))
        for k, x in enumerate(X):
            for i, _ in enumerate(self.classes):
                p = np.log(self.priors[i])
                for j, x_j in enumerate(x):
                    if abs(self.sigma[i][j]) < 1e-5:
                        continue
                    else:
                        p += np.log(self.gaussian(self.mu[i][j], self.sigma[i][j], x_j))
                proba[k, i] = np.exp(p)
        return proba

    def score(self, X, y):
        return np.mean(self.predict(X) == y)

    @staticmethod
    def gaussian(mu, sigma, x):
        return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (x - mu) * (x - mu) / (sigma * sigma))