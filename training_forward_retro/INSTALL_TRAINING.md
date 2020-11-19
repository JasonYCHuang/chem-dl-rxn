# INSTALL TRAINING


```bash
$ conda create --name retro-forward python=3.6
$ source activate retro-forward
```

```bash
$ conda install -c rdkit rdkit=2018.03.4.0
$ sudo apt-get install gcc libxrender1 libxext-dev
$ conda install pytorch torchvision cpuonly -c pytorch
```

```bash
$ git clone https://github.com/pytorch/fairseq
$ cd fairseq
$ pip install --editable ./
$ cd ..
```

```bash
$ git clone XXX
$ cd XXX
$ pip install -r requirements.txt
```
