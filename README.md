# tiny_md

> _tiny molecular dynamics_
>
> trabajo final de computación paralela 2020 (curso de posgrado - FaMAF, UNC) 

Proyecto de speed-up de un código de Dinámica Molecular con un potencial interatómico de Lennard-Jones. Inicialmente se dan las posiciones en un cristal FCC y velocidades aleatorias distribuidas uniformemente según la temperatura inicial. Las fuerzas se obtienen de un potencial de LJ (12-6) y la evolución temporal viene dada por el algoritmo Velocity Verlet. Se consideran condiciones periódicas de contorno para reproducir un sistema infinito. La temperatura es controlada a través de un reescaleo en las velocidades. Cada cierta cantidad de pasos de dinámica molecular se cambia la densidad del sistema y se reescalean las posiciones para obtener la ecuación de estado.

En el directorio `proyecto/` se encuentra el desarrollo del trabajo, para más información se puede ver `doc/informe.pdf`. 

En `proyecto0_cp2021/` está el proyecto inicial para la edición 2021 de Computación paralela.


### Requisitos

Para compilar es necesario tener instalado `gcc`, `OpenMP` y `OpenGL`.


### Compilación

Para correr la versión final del trabajo es necesario compilar utilizando `Makefile` en el directorio `src/`:
```bash
cd src/
make clean
make
```
donde `make clean` elimina los objetos compilados anteriormente y `make` compila dos ejecutables: `tiny_md` y `viz`, ambos ejecutables realizan la misma simulación pero el segundo posee una visualización en tiempo real.

> Nota:
>
> _Si se desean cambiar parámetros de entrada de la simulación, puede modificarse el archivo _`src/in.parameters`_ y recompilar o pasar los valores deseados como argumento a _`make`_; por ejemplo, _`make N=1372`_ cambia la cantidad de partículas que se simulan._


### Ejecución

Previo a la ejecución de `tiny_md` o `viz` puede configurarse OpenMP usando, por ejemplo, `OMP_PROC_BIND=true` y `OMP_NUM_THREADS=$i`, donde `$i` es la cantidad de hilos en los que se quiere correr.


### Contacto

Por errores, preguntas o sugerencias contactarse con:

+ Francisco Fernandez (<fernandezfrancisco2195@gmail.com>)
