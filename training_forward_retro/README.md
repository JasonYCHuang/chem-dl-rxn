# Preprocess

```
$ cd preprocess
$ bash prepare-rxn.sh
$ cd ..
```

```
# Preprocess/binarize the data
$ TEXT=preprocess/uspto
$ fairseq-preprocess --source-lang rct --target-lang prd \
    --trainpref $TEXT/train --validpref $TEXT/val --testpref $TEXT/test \
    --destdir data-bin/uspto.rct-prd \
    --workers 20
```

# Train

```
$ nohup fairseq-train \
    data-bin/uspto.rct-prd \
    --arch transformer_iwslt_de_en --share-decoder-input-output-embed \
    --optimizer adam --adam-betas '(0.9, 0.98)' --clip-norm 0.0 \
    --lr 5e-4 --lr-scheduler inverse_sqrt --warmup-updates 4000 \
    --dropout 0.3 --weight-decay 0.0001 \
    --criterion label_smoothed_cross_entropy --label-smoothing 0.1 \
    --max-tokens 4096 \
    --eval-bleu \
    --eval-bleu-args '{"beam": 5, "max_len_a": 1.2, "max_len_b": 10}' \
    --eval-bleu-detok moses \
    --eval-bleu-remove-bpe \
    --eval-bleu-print-samples \
    --best-checkpoint-metric bleu --maximize-best-checkpoint-metric &
```
