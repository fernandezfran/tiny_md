#include <GL/glut.h>    // OpenGL

#include "parameters.h"
#include "core.h"

// variables globales
static double Ekin, Epot, Temp, Pres;       // variables macroscopicas
static double Rho, cell_V, cell_L, tail, Etail, Ptail;
static double *rxyz, *vxyz, *fxyz;          // variables microscopicas
static double Rhob, sf, epotm, presm;
static int switcher = 0, frames = 0, mes;


// OpenGL specific drawing routines 

static double glL = cbrt((double)N / (Rhoi - 0.8));
static int win_id;
static int win_x = 600, win_y = 600;


static void pre_display ( void ){
	glViewport ( 0, 0, win_x, win_y );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
}


static void post_display ( void ){
	glutSwapBuffers ();
}


static void draw_atoms ( void ){
    
    int di;

    double dx;
    double dy;
    double dz;

    double reshape = 0.5 * ((glL - cell_L) / glL);

    glBegin( GL_POINTS );

        for (di = 0; di < 3*N; di+=3){ 
            // x, y, z entre 0 y 1
            dx = rxyz[di + 0] / glL + reshape;
            dy = rxyz[di + 1] / glL + reshape;
            dz = rxyz[di + 2] / glL + reshape;

            glColor3d(0.0, 1.0, 0.0);
            glVertex3d(dx, dy, dz);
        }

    glEnd ();

}


static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}


static void idle_func ( void ){
    
    if (switcher == 3){
        
        Rho = Rhoi;
        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail   = 16.0*M_PI*Rho*((2.0/3.0)*pow(rcut,-9) - pow(rcut,-3))/3.0;
        Etail  = tail*(double)N;
        Ptail  = tail*Rho;
        
        init_pos(rxyz, Rho);
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        switcher = 0;

    }else if (switcher == 2){ // imprimo propiedades en la terminal y cambio la densidad
        
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm/(double)mes,
                                   presm/(double)mes);

        Rhob   = Rho;
        Rho    = Rho - 0.1;


        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail   = 16.0*M_PI*Rho*((2.0/3.0)*pow(rcut,-9) - pow(rcut,-3))/3.0;
        Etail  = tail*(double)N;
        Ptail  = tail*Rho;
        
        sf = cbrt(Rhob/Rho);
        for (int k = 0; k < 3*N; k++){ // reescaleo posiciones a nueva densidad
            rxyz[k] *=sf;
        }
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);
        
        switcher = 0;
        if (fabs(Rho - (Rhoi - 0.9f)) < 1e-6){
            printf("\n");
            switcher = 3;
        }

    } else if (switcher == 1){ // loop de medición

        for (int i = frames; i < frames + tmes; i++){
            
            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, 
                            cell_V, cell_L);
            
            sf = sqrt(T0/Temp);
            for (int k = 0; k < 3*N; k++){ // reescaleo de velocidades
                vxyz[k] *=sf;
            }
        }
        
        Epot += Etail;
        Pres += Ptail;

        epotm += Epot;
        presm += Pres;
        mes++;

        frames += tmes;
        if (frames % trun == 0) switcher = 2;

    } else if (switcher == 0){ // loop de equilibración

        while ( frames % teq != 0){

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, 
                            cell_V, cell_L);

            sf = sqrt(T0/Temp);
            for (int k = 0; k < 3*N; k++){ // reescaleo de velocidades
                vxyz[k] *=sf;
            }

            frames++;
        }
        
        mes = 0;
        epotm = 0.0;
        presm = 0.0;

        switcher = 1;
    }
	glutSetWindow ( win_id );
	glutPostRedisplay ();
}


static void display_func ( void )
{
	pre_display ();
    draw_atoms ();
	post_display ();
}


static void open_glut_window ( void ){
    glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

    glutInitWindowPosition ( 0, 0 );
    glutInitWindowSize ( win_x, win_y );
    win_id = glutCreateWindow( "tiny molecular dynamics | visualization");

    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f);
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();

    pre_display ();
	
    // glutKeyboardFunc ( key_func );
	// glutMouseFunc ( mouse_func );
	// glutMotionFunc ( motion_func );
    
    glutReshapeFunc ( reshape_func );

    glutIdleFunc( idle_func );
	glutDisplayFunc ( display_func );
}


// viz main

int main( int argc, char ** argv){

    glutInit ( &argc, argv );

    rxyz = (double *)malloc(3*N*sizeof(double));
    vxyz = (double *)malloc(3*N*sizeof(double));
    fxyz = (double *)malloc(3*N*sizeof(double));

    // parametros iniciales para que los pueda usar (antes de modificar)
    // `idle_func`
    srand(SEED);
    Rho = Rhoi;
    Rhob   = Rho;
    cell_V = (double)N / Rho;
    cell_L = cbrt(cell_V);
    tail   = 16.0*M_PI*Rho*((2.0/3.0)*pow(rcut,-9) - pow(rcut,-3))/3.0;
    Etail  = tail*(double)N;
    Ptail  = tail*Rho;
    
    init_pos(rxyz, Rho);
    init_vel(vxyz, &Temp, &Ekin);
    forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);
    //
    //

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", teq);
    printf("# Pasos de medición:         %d\n", trun - teq);
    printf("# (mediciones cada %d pasos)\n", tmes);
   
    open_glut_window ();

    glutMainLoop();

    exit ( 0 );
}
