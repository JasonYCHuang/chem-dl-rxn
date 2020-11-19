#!/usr/bin/env bash
#
# Adapted from https://github.com/facebookresearch/MIXER/blob/master/prepareData.sh

MOSES_GIT=https://github.com/moses-smt/mosesdecoder.git
if [ -d 'mosesdecoder' ]; then
    echo "Mosesdecoder has already been cloned."
else
    echo 'Cloning Moses github repository (for tokenization scripts)...'
    git clone $MOSES_GIT
fi

SUBWORD_GIT=https://github.com/rsennrich/subword-nmt.git
if [ -d 'subword-nmt' ]; then
    echo "Subword-nmt has already been cloned."
else
    echo 'Cloning Subword NMT repository (for BPE pre-processing)...'
    git clone $SUBWORD_GIT
fi



BPEROOT=subword-nmt/subword_nmt
BPE_TOKENS=600

GZ=download_retro
prep=uspto
tmp=$prep/tmp
orig=orig

src=prd
tgt=rct
# condition=MIT_mixed

if [ -d 'orig' ]; then
    echo "Data has already been downloaded."
else
    mkdir -p $orig $tmp $prep
    echo "Downloading data from ${URL}..."
    cd $orig
    cp '../../backup/download_retro' .

    if [ -f $GZ ]; then
        echo "Data successfully downloaded."
    else
        echo "Data not successfully downloaded."
        exit
    fi
    tar zxvf $GZ
    cd ..
fi

TRAIN=$tmp/train.prd-rct
BPE_CODE=$prep/code
rm -f $TRAIN
for l in $src $tgt; do
    cat $orig/$prep/$condition/train.$l >> $TRAIN
done

echo "learn_bpe.py on ${TRAIN}..."
python $BPEROOT/learn_bpe.py -s $BPE_TOKENS < $TRAIN > $BPE_CODE

for L in $src $tgt; do
    for f in train.$L val.$L test.$L; do
        echo "apply_bpe.py to ${f}..."
        python $BPEROOT/apply_bpe.py -c $BPE_CODE < $orig/$prep/$condition/$f > $prep/$f
    done
done

echo "BPE_TOKENS=${BPE_TOKENS}"

echo "DONE!..."


