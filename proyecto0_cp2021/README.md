# tiny_md

> _tiny molecular dynamics_

Proyecto inicial para realizar speed-up de un código de Dinámica Molecular con un potencial interatómico de Lennard-Jones. Inicialmente se dan las posiciones en un cristal FCC y velocidades aleatorias distribuidas uniformemente según la temperatura inicial. Las fuerzas se obtienen de un potencial de LJ (12-6) y la evolución temporal viene dada por el algoritmo Velocity Verlet. Se consideran condiciones periódicas de contorno para reproducir un sistema infinito. La temperatura es controlada a través de un reescaleo en las velocidades. Cada cierta cantidad de pasos de dinámica molecular se cambia la densidad del sistema y se reescalean las posiciones para obtener la ecuación de estado. Para más información se puede ver `informe.pdf`.


### Requisitos

Para compilar es necesario tener instalado `gcc` y `OpenGL`.


### Compilación

Para compilar se utiliza `Makefile`:
```bash
make clean
make
```
donde `make clean` elimina los objetos compilados anteriormente y `make` compila dos ejecutables: `tiny_md` y `viz`, ambos realizan la misma simulación pero el segundo posee una visualización en tiempo real.

> Nota:
>
> _Si se desean cambiar parámetros de entrada de la simulación, puede modificarse el archivo _`in.parameters`_ y recompilar o pasar los valores deseados como argumento a _`make`_; por ejemplo, _`make N=1372`_ cambia la cantidad de partículas que se simulan._


### Contacto

Por errores, preguntas o sugerencias contactarse con:

+ Francisco Fernandez (<fernandezfrancisco2195@gmail.com>)
