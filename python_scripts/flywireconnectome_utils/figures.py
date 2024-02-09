from matplotlib import pyplot as plt


def savefig(path, dpi=300):
    if path.endswith(".png"):
        plt.savefig(path, dpi=dpi)
        plt.savefig(f"{path[:-4]}.pdf")
