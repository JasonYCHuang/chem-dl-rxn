def translate_top_k(model, inputs, k=5):
    in_bins = []
    for inp in inputs:
        toks = model.tokenize(inp)
        bpe = model.apply_bpe(toks)
        in_bin = model.binarize(bpe)
        in_bins.append(in_bin)

    out_bins = model.generate(in_bins, beam=k)
    results = []
    for ob in out_bins:
        outs = []
        for o in ob:
            bpe = model.string(o['tokens'])
            toks = model.remove_bpe(bpe)
            out = model.detokenize(toks)
            outs.append(out)
        results.append(list(set(outs))[:k])
    return results
