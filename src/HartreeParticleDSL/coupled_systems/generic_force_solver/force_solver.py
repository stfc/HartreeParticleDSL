from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
from HartreeParticleDSL.HartreeParticleDSLExceptions import CoupledSystemError

class force_solver(base_coupler):
    def __init__(self):
        pass

    def get_includes(self):
        return []

    def get_includes_header(self):
        return []

    def MPI_allowed(self, MPI_status):
        if MPI_status:
            raise CoupledSystemError(f"{type(self)} doesn't support MPI.")

