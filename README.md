# LoRP-FDM-DCD
## A mixed-integer linear programming model and a matheuristic for the location-or-routing problem with fixed destination multi-depot and distance constrained drones
This code was written in the paper ``An optimization-based approach for an integral forest fires monitoring system with multiple technologies and surveillance drones: A case study in Southern Chile'' by Rodrigo De la Fuente and Maichel M. Aguayo and Carlos Contreras-Bolton.

The instances are in the directory `instances/`, both sets of instances, random instances and based on a real-life case study in Chile.

## Instructions to run

To build the binary type:

```
make
```
To run the program:

```
./run.sh algorithm instance [instance_i] [instance_j] [time] [memory]

```
where **algorithm** is:
```
M3-MTZ: proposed model with MTZ formulation
M3-GG: proposed model with GG formulation
MH: matheuristic
M1: variant M1 of the matheuristic
M2: variant M2 of the matheuristic
```

where **instance** is:
```
A: random instances
B: instances based on a real-life case study in Chile
```

where **instance_i** and **instance_j** are optional and mean the range of the instance, for A [1, 10] and B [1, 2].
where **time** is the time limit in seconds and **memory** is the memory limit in Megabytes.

Give an example
```
/run.sh MH A 1 10 300 5000
```

## License
The program is distributed under the GGNU Affero General Public License v3.0
See the `LICENSE` file.
