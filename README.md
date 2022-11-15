
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1667380285/logos/pygomodo_logo_duolry.png"  width="400" align='center' style="margin-top:0px;margin-left:150px"/>

<br>
pyGOMoDo is a Python library to perform homology modeling and docking specifically designed for human GPCRs. pyGOMoDo is a python wrap-up of the updated functionalities of GOMoDo web-service (https://gomodo.grs.kfa-juelich.de). It was developed having in mind the ist usage through jupyter notebook, where users can create their own protocols of modeling and docking of GPCRs.

<br><br>

# Quickstart
We recomend to use pyGOMoDo with Docker([here](https://hub.docker.com/r/rpribeiro/pygomodo)). The container deploys the pyGOMoDo jupyter environment.

### Running the tutorials

Linux :
```bash
$ docker run -d -t -p 8886:8885 --name pygomodo_tutorials rpribeiro/pygomodo
```
Then open the the following link on the browser:
```
localhost:8886
```
<br>

### Running production environment (reading and writing local directories)
Linux :
```bash
$ docker run -d -t -p 8886:8885 -v $PWD:$PWD:Z -w $PWD --name pygomodo_production rpribeiro/pygomodo
```
Then open the the following link on the browser:
```
localhost:8886
```
<br>

## Documentation

All documentation [here](https://pygomodo.readthedocs.io).

<br>

# Developed by
<div style="padding-bottom:50px">
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1637335206/logos/Logo_des_Forschungszentrums_J_C3_BClich_seit_2018_hcliq4.svg"  width="200" align='left' style="margin-top:10px"/>
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1667384566/logos/empty_sxac7h.png"  width="50" align='left' style="margin-top:10px"/>
<img src="https://res.cloudinary.com/djz27k5hg/image/upload/v1657885120/logos/univr_logo_rspn8o.jpg"  width="200" style="margin-top:10px; margin-left:5000px"/>


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
@article{XXX,
    title={XXX},
    author={XXX},
    publisher={XXX},
    note={\url{XXX}},
    year={XXX}
}
```