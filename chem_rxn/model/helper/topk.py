from rdkit import Chem

def translate_top_k(model, inputs, k=5):
    inputs = [Chem.MolToSmiles(Chem.MolFromSmiles(inp)) for inp in inputs]

    in_bins = []
    for inp in inputs:
        toks = model.tokenize(inp)
        bpe = model.apply_bpe(toks)
        in_bin = model.binarize(bpe)
        in_bins.append(in_bin)

    out_bins = model.generate(in_bins, beam=k)
    results = []
    for idx, ob in enumerate(out_bins):
        outs = []
        for o in ob:
            bpe = model.string(o['tokens'])
            toks = model.remove_bpe(bpe)
            out = model.detokenize(toks)
            if out not in outs:
                outs.append(out)
        results.append(outs[:k])
    return results
