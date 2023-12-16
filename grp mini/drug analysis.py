from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS

def detect_isomorphism(mol1, mol2):
    # Convert molecules to Mol objects
    mol1 = Chem.MolFromSmiles(mol1)
    mol2 = Chem.MolFromSmiles(mol2)
    
    # Generate MCS (Maximum Common Substructure)
    mcs = rdFMCS.FindMCS([mol1, mol2])
    query = Chem.MolFromSmarts(mcs.smartsString)
    
    # Check if molecules contain the MCS substructure
    matches1 = mol1.GetSubstructMatches(query)
    matches2 = mol2.GetSubstructMatches(query)
    
    # Highlight and display the connected substructure
    if len(matches1) > 0 and len(matches2) > 0:
        substructure1 = Chem.MolFromSmarts(mcs.smartsString)
        substructure2 = Chem.MolFromSmarts(mcs.smartsString)
        
        # Highlight the substructure in molecule 1
        for match in matches1:
            for idx in match:
                substructure1.GetAtomWithIdx(idx).SetProp("atomNote", "matched")
        
        # Highlight the substructure in molecule 2
        for match in matches2:
            for idx in match:
                substructure2.GetAtomWithIdx(idx).SetProp("atomNote", "matched")
        
        # Display the molecules with the connected substructure highlighted
        img = Draw.MolsToGridImage([mol1, mol2, substructure1, substructure2], molsPerRow=2,
                                   subImgSize=(300, 300), legends=["Molecule 1", "Molecule 2", "Substructure 1", "Substructure 2"])
        img.show()
    else:
        print("No isomorphic substructure found.")

# Example usage
mol1 = "CC(=O)OC1=CC=CC=C1C(=O)O"
mol2 = "CC(=O)OC1=CC=C(C=C1)C(=O)O"
detect_isomorphism(mol1, mol2)
