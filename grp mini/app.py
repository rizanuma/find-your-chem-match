from rdkit import Chem
from rdkit.Chem import Draw
from difflib import SequenceMatcher
from flask import Flask, request, render_template
import base64

app = Flask(__name__)

def calculate_similarity(smiles1, smiles2):
    seq_matcher = SequenceMatcher(None, smiles1, smiles2)
    return seq_matcher.ratio() * 100


def analyze_drug(input_smiles):
    stored_drugs = ["CC(=O)OC1=CC=CC=C1C(=O)O",
                     "CC(=O)NC1=CC=C(C=C1)O", 
                     "CC(=O)OC1=C(C=CC=C1O)C(=O)O",
                     "F/C=C/F or F\C=C\F ",
                     "C[C@@H](C(=O)O)N",
                     ]

    best_match = None
    best_similarity = 0

    for drug in stored_drugs:
        if Chem.MolToSmiles(Chem.MolFromSmiles(input_smiles)) == Chem.MolToSmiles(Chem.MolFromSmiles(drug)):
            best_match = drug
            best_similarity = 100
            break
        else:
            similarity = calculate_similarity(input_smiles, drug)
            if similarity > best_similarity:
                best_match = drug
                best_similarity = similarity

    return best_match, best_similarity

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        input_smiles = request.form['smiles']
        best_match, best_similarity = analyze_drug(input_smiles)
        mol = Chem.MolFromSmiles(best_match)
        image = Draw.MolToImage(mol)
        image_data = image_to_base64(image)

        html = f'''
        <h1>Drug Analysis</h1>
        <form method="POST" action="/">
            <input type="text" name="smiles" placeholder="Enter Drug SMILES" value="{input_smiles}">
            <input type="submit" value="Analyze">
        </form>
        
        <h2>Results:</h2>
        <ul>
            <li>Drug: {best_match}</li>
            <li>Similarity: {best_similarity}%</li>
            <li>Structure:</li>
            <img src="data:image/png;base64,{image_data}" alt="Drug Structure">
        </ul>
        '''

        return html
    else:
        html = '''
        <h1>Drug Analysis</h1>
        <form method="POST" action="/">
            <input type="text" name="smiles" placeholder="Enter Drug SMILES">
            <input type="submit" value="Analyze">
        </form>
        '''

        return html

def image_to_base64(image):
    import io
    import base64

    buffered = io.BytesIO()
    image.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return img_str

if __name__ == '__main__':
    app.run(debug=True)
