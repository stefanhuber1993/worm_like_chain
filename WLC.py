import numpy as np
import matplotlib.pyplot as plt


class ChainBundle():
    def __init__(self):
        self.data = None

    def compute(self, n_chains, persistence_length, contour_length, n_segments=100):
        """
        :param n_chains: How many chains to include in the bundle
        :param persistence_length: Average correlation between chain segments gets down to 1/e after persistence_length
        :param contour_length: Length of the Chains to be simulated
        :param n_segments: How smooth chain should be approximated by linear segments.
                           Should be at least 50 for meaningful results
        :writes in self.data: np.array of shape (n_chains, n_segments, 2) (last entry contains x and y coordinate)
        """
        segment_length = contour_length / n_segments
        persistence_length = float(persistence_length)

        coordinates = np.zeros((n_chains, n_segments, 2), dtype=np.float64)  # Initialise numpy array
        coordinates[:, 0, :] = 0  # Chains start at (0,0)
        coordinates[:, 1, 0] = segment_length  # Chain starts directly into x direction
        angle = [0]

        for segment in range(n_segments - 2):
            random_dir = np.random.choice([-1, 1], n_chains).astype(np.float64)
            angle += random_dir * np.arccos(np.exp(- segment_length / persistence_length / 1.0))
            coordinates[:, segment + 2, 0] = coordinates[:, segment + 1, 0] + segment_length * np.cos(angle)
            coordinates[:, segment + 2, 1] = coordinates[:, segment + 1, 1] + segment_length * np.sin(angle)

        self.data = coordinates

    def compute_persistence_length(self, verbose=False):
        """
        verbose == True:
        :return: persistence length, x_values, exponential_curve
        verbose == False:
        :return: persistence length
        """
        from scipy.optimize import curve_fit
        """Similar to
        https://pythonhosted.org/MDAnalysis/_modules/MDAnalysis/analysis/polymer.html"""
        n = self.data.shape[1]
        results_total = np.zeros((self.data.shape[0], n - 1))
        for a, chain in enumerate(self.data):
            chain = chain
            vecs = chain[1:] - chain[:-1]
            vecs_norm = vecs / np.sqrt((vecs * vecs).sum(axis=1))[:, None]
            inner_pr = np.inner(vecs_norm, vecs_norm)
            results = np.zeros(n - 1)
            for i in range(n - 1):
                results[:(n - 1) - i] += inner_pr[i, i:]
            norm = np.linspace(n - 1, 1, n - 1)
            results = results / norm
            results_total[a] = results
        results_curve = results_total.mean(0)

        def expon(lamb, x):
            return np.exp(-lamb * x)

        l = np.sqrt(np.sum(vecs[0] ** 2))
        x = np.arange(n - 1) * l
        lamb, dlamb = curve_fit(expon, x, results_curve, p0=0.1)
        if verbose:
            return 1 / lamb[0], x, results_curve
        else:
            return 1 / lamb[0]

    def plot(self, ax=None, **kwargs):
        """
        :param kwargs: e.g. color='blue'
        :return: ax
        """
        if ax is None:
            ax = plt.gca()
        label = kwargs.pop('label','')
        if self.data is None:
            raise "call .compute first to get a WLC bundle of chains"
        else:
            for i, chain in enumerate(self.data):
                ax.plot(*chain.T, label="%s"%label if i == 0 else '', **kwargs)
        return ax
