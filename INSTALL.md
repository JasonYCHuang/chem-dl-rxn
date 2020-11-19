# INSTALL


```bash
$ conda create --name retro-forward python=3.6
$ source activate retro-forward
```

```bash
$ conda install -c rdkit rdkit=2020.03.3.0
$ sudo apt-get install gcc libxrender1 libxext-dev libxml2-dev libxslt1-dev zlib1g-dev g++ unzip
$ conda install pytorch torchvision cpuonly -c pytorch
$ conda install -c anaconda flask
```

```bash
$ git clone https://github.com/pytorch/fairseq
$ cd fairseq
$ pip install --editable ./
$ cd ..
```


# Get `chem-rxn-app`


```bash
$ git clone git@github.com:JasonYCHuang/chem-rxn-app.git
$ cd chem-rxn-app
$ pip install -r requirements.txt
```


# Get the Model

```bash
$ wget https://bwsyncandshare.kit.edu/s/NGG6oqCJkdGnHx4/download
$ unzip download
$ mv nlp-rxn-main chem_rxn
$ rm download
```

# Run

```bash
# run on the production server
$ gunicorn -w 4 -b 0.0.0.0:7007 server:app --daemon
```


```bash
# for local development only
$ export FLASK_APP=chem_rxn && export FLASK_ENV=development && flask run --host=0.0.0.0 --port=7007
```
