# pyEXP-examples/Tutorials/EXP

Scripts and notebooks that give examples of running the `exp` N-body
code and analyzing it briefly using pyEXP

These simulations use your installed `exp` installation. If you have
not build `exp` but only `pyEXP`, you can try the Docker image.  `exp`
and the necessary standalone support routines will be automatically
available in your Docker container.

The `exp` configuration files, `config.yml` included in each example, assume
the default library location for the Docker container:
```
/usr/local/EXP/lib/user
```
You will need to adjust this for a local installation of `exp`, e.g.:
```
/home/user/lib/user
```

These simulations use a very small number of particles and are
intended to be used as a starting point for learning how to use
`exp`. Nonetheless, they will take will take many minutes (up to an
hour) to run on a laptop.

## Organization

| Directory    | Contents |
| ---          | ---      |
| DiskHalo     | Generate disk and halo initial conditions and run a galaxy simulation |
| Cube         | Run a Jeans unstable, initially homogeneous distribution in a periodic cube |
| Zang         | Run an unstable, doubly tapered Mestel disk to produce the classic m=2 spiral |
