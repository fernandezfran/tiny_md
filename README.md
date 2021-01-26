Proyecto de speed-up de un código de Dinámica Molecualar NVT (reescaleo de velocidades) para un fluido de Lennard-Jones como parte del curso "Computación Paralela 2020" (https://cs.famaf.unc.edu.ar/~nicolasw/Docencia/CP/2020/index.html)

* 0-proyecto_inical/

    código inicial con Makefile y parámetros de entrada (flags gcc = https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)

* 1-optimizacion_secuencial

    1. 01-flags-aos/

        exploración de flags sin realizar ningún cambio en el proyecto inicial

    2. 02-flags-soa/

        cambio el orden de almacenamiento de los datos en los vectores, en vez de guardar x, y, z (para pos, vel y frc) consecutivos los separo en vectores de largo N (es decir, es equivalente a tener una matriz 3xN, donde en la fila 1 está x, en la 2 y, en la 3 z).

        (es decir, antes de tenía AoS, ahora SoA).

    3. 03-mixed-precision/

        veo si en algún lugar puedo mezclar precisiones, usar float en vez de double

* 2-vectorizacion/
 
    1. autovectorizacion/

        ¿qué loops autovectoriza gcc de acuerdo a si tengo SoA o AoS? ¿ayuda agregar loops naiv?

    2. sse/

        vectorización del cálculo de fuerzas con sse (empezado, falta).
