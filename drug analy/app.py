from flask import Flask, request,render_template
import base64
from io import BytesIO
from rdkit.Chem.Draw import MolToImage
from flask import render_template_string
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

app = Flask(__name__)

def draw_molecule(mol):
    image = MolToImage(mol)
    buffered = BytesIO()
    image.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return img_str


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        smiles1 = request.form["smiles1"]
        smiles2 = request.form["smiles2"]
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        img_mol1 = draw_molecule(mol1)
        img_mol2 = draw_molecule(mol2)

        return render_template_string(
            '''
            <h1>Molecule Similarity Comparison</h1>
            <h2>Molecule 1:</h2>
            <img src="data:image/png;base64, {{ img_mol1 }}" />
            <h2>Molecule 2:</h2>
            <img src="data:image/png;base64, {{ img_mol2 }}" />
            ''',
            img_mol1=img_mol1,
            img_mol2=img_mol2
        )

    return render_template("index.html")

if __name__ == '__main__':
    app.run()
