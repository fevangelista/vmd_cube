import pyscf, pyscf.tools
import numpy as np

def main():
    mol = pyscf.gto.M(atom="""
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """, basis='sto-3g',symmetry=True)

    mf = pyscf.scf.RHF(mol)
    mf.kernel()

    if (mol.symmetry == False):
        for i in range(mol.nao):
            pyscf.tools.cubegen.orbital(mol, f'Psi_a_{i}_{i}-A.cube', mf.mo_coeff[:,i], resolution=0.2, margin=3.0)
    else:
        irrep = mf.mo_coeff.orbsym
        irrep_indices = irrep.copy()

        for i in mol.irrep_id:
            irrep_indices[np.argwhere(irrep==i).flatten()] = np.arange(1, 1+np.sum(irrep==i))

        irrep_dict = {mol.irrep_id[i]:mol.irrep_name[i] for i in range(len(mol.irrep_id))}

        for i in range(mol.nao):
            pyscf.tools.cubegen.orbital(mol, f'Psi_a_{i+1}_{irrep_indices[i]}-{irrep_dict[mf.mo_coeff.orbsym[i]]}.cube', mf.mo_coeff[:,i], resolution=0.2, margin=3.0)

if __name__ == '__main__':
    main()