import pyscf, pyscf.tools
import numpy as np

def main():
    mol = pyscf.gto.M(atom="""
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """, basis='sto-3g',symmetry='c2v')

    mf = pyscf.scf.RHF(mol)
    mf.kernel()

    irrep = mf.mo_coeff.orbsym
    irrep_indices = irrep.copy()

    irrep_indices[np.argwhere(irrep==0).flatten()] = np.arange(1, 1+np.sum(irrep==0))
    irrep_indices[np.argwhere(irrep==1).flatten()] = np.arange(1, 1+np.sum(irrep==1))

    irrep_dict = {mol.irrep_id[i]:mol.irrep_name[i] for i in range(len(mol.irrep_id))}

    for i in range(mol.nao):
        pyscf.tools.cubegen.orbital(mol, f'Psi_a_{i}_{irrep_indices[i]}-{irrep_dict[mf.mo_coeff.orbsym[i]]}.cube', mf.mo_coeff[:,i])

if __name__ == '__main__':
    main()