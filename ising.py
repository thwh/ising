import click
from dataclasses import dataclass
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import numpy as np
from tqdm import tqdm


@dataclass
class IsingLattice:
    size: int
    T: float
    system: np.array


def BuildSystem(initial_state, size, temperature):
    """Build the system

    Build either a randomly distributed system or a homogeneous system (for
    watching the deterioration of magnetization

    Parameters
    ----------
    initial_state : str
        Initial state of the lattice.  currently only random 'r' initial
        state, or 'u' uniformly magnetized, is supported
    size : int
        Size of the lattice.
    temperature : float
        Global temperature of the system.

    Return
    ------
    IsingLattice
        dataclass instance
    """
    sq_size = (size, size)

    if initial_state == 'r':
        a = np.random.choice([-1, 1], sq_size)
    elif initial_state == 'u':
        a = np.ones(sq_size)
    else:
        raise ValueError(
            "Initial State must be 'r', random, or 'u', uniform"
        )
    return IsingLattice(size, temperature, a)


def BoundaryCondition(i, lattice):
    """Apply periodic boundary condition

    Check if a lattice site coordinate falls out of bounds. If it does,
    apply periodic boundary condition

    Parameters
    ----------
    i : int
        lattice site coordinate
    lattice : np.array
        Nested lists representing the 2D lattice.

    Return
    ------
    int
        corrected lattice site coordinate
    """
    return i % lattice.size


def CalculateEnergy(lattice, N, M):
    """Calculate the energy of spin interaction at a given lattice site
    i.e. the interaction of a Spin at lattice site n,m with its 4 neighbors

    - S_n,m*(S_n+1,m + Sn-1,m + S_n,m-1, + S_n,m+1)

    Parameters
    ----------
    lattice : np.array
        Nested lists representing the 2D lattice.
    N : int
        lattice site coordinate
    M : int
        lattice site coordinate

    Return
    ------
    float
        energy of the site
    """

    neighbor_energy = 0
    for n in (N-1, N+1):
        neighbor_energy += lattice.system[BoundaryCondition(n, lattice), M]
    for m in (M-1, M+1):
        neighbor_energy += lattice.system[N, BoundaryCondition(m, lattice)]
    return -2 * lattice.system[N, M] * neighbor_energy


def UpdateLattice(lattice):
    # Randomly select a site on the lattice
    N, M = np.random.randint(0, lattice.size, 2)

    # Calculate energy delta of a flipped spin
    ΔE = -1*CalculateEnergy(lattice, N, M)

    # if energetic change is dissipative
    if ΔE <= 0.0:
        lattice.system[N, M] *= -1

    # Or, if e^(-ΔE*β) > some [0, 1]
    #
    # Here β is the thermodynamic beta, which is 1/kT, where k is the
    # Boltzmann constant.
    #
    # That constant is thrown out when not dealing with molecular
    # movements. So we're left with e^(-ΔE * 1/T):
    elif np.exp(-ΔE * 1/lattice.T) > np.random.rand():
        lattice.system[N, M] *= -1

    return lattice


def InternalEnergy(lattice):
    """Calculate the internal energy of the system.

    Parameters
    ----------
    lattice : np.array
        Nested lists representing the 2D lattice.

    Return
    ------
    float, float
        energy and first exponent
    """

    e = 0
    E = 0
    E_2 = 0

    for i in range(lattice.size):
        for j in range(lattice.size):
            e = CalculateEnergy(lattice, i, j)
            E += e
            E_2 += e**2

    U = (1/lattice.size**2)*E
    U_2 = (1/lattice.size**2)*E_2

    return U, U_2


def HeatCapacity(lattice):
    U, U_2 = InternalEnergy(lattice)
    return U_2 - U**2


def Magnetization(lattice):
    """Find the overall magnetization of the system
    """
    return np.abs(np.sum(lattice.system)/lattice.size**2)


def Run(lattice, epochs, video=True):
    """Run the simulation
    """

    FFMpegWriter = manimation.writers['ffmpeg']
    writer = FFMpegWriter(fps=10)

    fig = plt.figure()

    with writer.saving(fig, "ising.mp4", dpi=100):
        for epoch in tqdm(range(epochs)):

            lattice = UpdateLattice(lattice)

            if video and epoch % (epochs//75) == 0:
                img = plt.imshow(
                    lattice.system, interpolation='nearest', cmap='jet'
                )
                writer.grab_frame()
                img.remove()

    plt.close('all')


@click.command()
@click.option(
    '--temperature', '-t',
    default=0.5,
    show_default=True,
    help='temperature of the system'
)
@click.option(
    '--initial-state', '-i',
    default='r',
    type=click.Choice(['r', 'u'], case_sensitive=False),
    show_default=True,
    help='(R)andom or (U)niform initial state of the system'
)
@click.option(
    '--size', '-s',
    default=100,
    show_default=True,
    help='Number of sites, M, in the MxM lattice'
)
@click.option(
    '--epochs', '-e',
    default=1_000_000,
    type=int,
    show_default=True,
    help='Number of iterations to run the simulation for'
)
@click.option(
    '--video', '-v',
    is_flag=True,
    default=True,
    help='Record a video of the simulation progression'
)
def main(temperature, initial_state, size, epochs, video):
    lattice = BuildSystem(
        initial_state=initial_state, size=size, temperature=temperature
    )
    Run(lattice, epochs, video)

    print(f"{'Net Magnetization [%]:':.<25}{Magnetization(lattice):.2f}")
    print(f"{'Heat Capacity [AU]:':.<25}{HeatCapacity(lattice):.2f}")


if __name__ == "__main__":
    plt.ion()
    main()
