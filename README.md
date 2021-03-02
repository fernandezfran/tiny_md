# tiny_md

> _tiny molecular dynamics_
>
> > trabajo final de computación paralela 2020 (curso de posgrado - FaMAF, UNC) 

Proyecto de speed-up de un código de Dinámica Molecular con un potencial interatómico de Lennard-Jones. En el directorio `proyecto/` se encuentra el desarrollo del trabajo, para más información ver `doc/informe.pdf`.


### Requisitos

Para compilar es necesario tener instalado `gcc`, `OpenMP` y `OpenGL`.


### Compilación

Para correr la versión final del trabajo es necesario compilar utilizando `Makefile` en el directorio `src/`:
```bash
cd src/
make clean
make
```
donde `make clean` elimina los objetos compilados anteriormente y `make` compila dos ejecutables: `tiny_md` y `tiny_md_viz`, el segundo posee una visualización en tiempo real.

> Nota:
>
> > _Si se desean cambiar parámetros de entrada de la simulación, puede modificarse el archivo _`src/in.parameters`_ y recompilar o pasar los valores deseados como argumento a _`make`_; por ejemplo, _`make N=1372`_ cambia la cantidad de partículas que se simulan._


### Ejecución

Previo a la ejecución de `tiny_md` o `tiny_md_viz` puede configurarse OpenMP usando, por ejemplo, `OMP_PROC_BIND=true` y `OMP_NUM_THREADS=$i`, donde `$i` es la cantidad de hilos en los que se quiere correr.
