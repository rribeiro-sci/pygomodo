<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1667380285/logos/pygomodo_logo_duolry.png"  width="400" align='center' style="margin-top:0px;margin-left:150px"/>

<br><br>

[![DOI](https://zenodo.org/badge/497925678.svg)](https://zenodo.org/badge/latestdoi/497925678)
[![Documentation Status](https://readthedocs.org/projects/pygomodo/badge/?version=latest)](https://pygomodo.readthedocs.io/en/latest/?badge=latest)
![GitHub](https://img.shields.io/github/license/rribeiro-sci/pygomodo)

![GitHub release (latest by date)](https://img.shields.io/github/v/release/rribeiro-sci/pygomodo)
![Docker Image Size (latest by date)](https://img.shields.io/docker/image-size/rpribeiro/pygomodo)
![Docker Pulls](https://img.shields.io/docker/pulls/rpribeiro/pygomodo)
![Docker Stars](https://img.shields.io/docker/stars/rpribeiro/pygomodo)

# Description
pyGOMoDo is a Python library to perform homology modeling and docking specifically designed for human GPCRs. pyGOMoDo is a python wrap-up of the updated functionalities of GOMoDo web-service (https://gomodo.grs.kfa-juelich.de). It was developed having in mind the ist usage through jupyter notebook, where users can create their own protocols of modeling and docking of GPCRs.

<br><br>

# Quickstart
We recomend to use pyGOMoDo with Docker([here](https://hub.docker.com/r/rpribeiro/pygomodo)). The container deploys the pyGOMoDo jupyter environment.


### Running pyGOMoDo docker container (reading and writing local directories)

```bash
docker run -d -t -p 8886:8885 -v $PWD:$PWD:Z --name pygomodo_production rpribeiro/pygomodo
```
Then open the the following link on the browser:
```
localhost:8886
```
<br>

### How to install docker

Instructinos to install docker [here](https://github.com/rribeiro-sci/pygomodo/wiki/How-to-install-Docker).


### Running pyGOMoDo locally

The instructions to install pyGOMoDO [here](https://github.com/rribeiro-sci/pygomodo/wiki/How-to-install-Docker).

<br>


## Documentation
Documentation can be found online on [ReadTheDocs](https://pygomodo.readthedocs.io). 

<br>

# Developed on behalf of
<div style="padding-bottom:50px">
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1637335206/logos/Logo_des_Forschungszentrums_J_C3_BClich_seit_2018_hcliq4.svg"  width="200" align='left' style="margin-top:10px"/>
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1667384566/logos/empty_sxac7h.png"  width="50" align='left' style="margin-top:10px"/>
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1657885120/logos/univr_logo_rspn8o.jpg"  width="200" style="margin-top:10px; margin-left:1px"/>


<br><br><br><br><br>

# Funded by
EU Human Brain Project (SGA1 and SGA2): This open source software was developed in part in the Human Brain Project funded from the European Union's Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No 720270 and No. 78907 (Human Brain Project SGA1 and SGA2).
<div style="padding-bottom:50px">
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1637657234/logos/HBP_horizontal_logo_qtcyzn.png" width="300" align='left' style="margin-left:50px">
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1642677502/logos/COFUNDED_EU_j2ktlp.jpg" width="300" align='left' style="margin-left:50px">
</div> 

<br><br><br><br>


# Cite Us
If you use or adapt the pyGOMoDo for your own research projects please cite us.

```
@article{10.1093/bioinformatics/btad294,
    author = {Ribeiro, Rui P and Giorgetti, A},
    title = "{pyGOMoDo: GPCRs modeling and docking with python}",
    journal = {Bioinformatics},
    year = {2023},
    doi = {10.1093/bioinformatics/btad294},
    url = {https://doi.org/10.1093/bioinformatics/btad294},
}
```
