//an improved, yet simple 3d gravity simulator

#include<stdio.h>
#include<GL/glut.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define WIDTH	 	750
#define HEIGHT 		750
#define NUMTRACE 	50
#define TIMESTEP 	20

#define TEXT_RED	"\033[31m"
#define TEXT_YELLOW "\033[33m"
#define TEXT_GREEN	"\033[32m"
#define TEXT_WHITE 	"\033[37m"


//planet structure
typedef struct {
	double pos[ 3 ];
	double vel[ 3 ];
	double acc[ 3 ];
	double rad;
	double mass;
	double netforce[ 3 ];
	double trace[ NUMTRACE ][ 3 ];
	float col[ 3 ];
	_Bool collided;
	_Bool tracing;
} obj_t;

//universe structure
typedef struct {
	obj_t *bodies;
	int nBodies;
	int time;
	float zoom, alpha[ 1000 ];
	double azi, alt, zen;
	_Bool tracing;
} uni_t;
uni_t *sim; //this has to be global for GLUT

//program functions
double randDouble( double min, double max );
void zero( int nBodies );
void save( void );
void load( void );
void simulate( void );
void draw( void );
void keys( unsigned char key, int x, int y );
void timer( );

double randDouble( double min, double max ) {
    double random = (( double ) rand( ) ) / ( double ) RAND_MAX;
    double range = max - min;
    return( random * range ) + min;
}

//kinda like the big bang
void zero( int nBodies ) {
	sim = ( uni_t* )malloc( sizeof( uni_t ) );
	sim->nBodies = nBodies;
	sim->bodies = ( obj_t* )malloc( sizeof( obj_t ) * sim->nBodies );
	sim->time = 0;
	sim->zoom = 1.0;
	for( int i = 0; i < NUMTRACE; i ++ ) {
		sim->alpha[ i ] = ( float )( i + 1 ) / NUMTRACE;
	}
	sim->azi = 0.0;
	sim->alt = 0.0;
	sim->zen = 0.0;

	sim->tracing = 0;

	obj_t *this = sim->bodies;
	for( int i = 0; i < sim->nBodies; i ++ ) {
		this->rad = randDouble( 0.001, 0.01 );
		this->mass = 4.19 * this->rad;
		for( int j = 0; j < 3; j ++ ) {
			this->pos[ j ] = randDouble( 0.0, 2.0 ) - 1.0;
			this->col[ j ] = randDouble( 0.0, 1.0 );
		}
		for( int j = 0; j < NUMTRACE; j ++ ) {
			for( int k = 0; k < 3; k ++ ) {
				this->trace[ j ][ k ] = this->pos[ k ];
			}
		}
		this ++; //increment body pointer
	}
	printf( "Initialised simulation with %i bodies.\n", nBodies );
	return;
}

void keys( unsigned char key, int x, int y ) {
	int win;
	switch( key ) {
		case 27:
			win = glutGetWindow( );
			glutDestroyWindow( win );
			exit( EXIT_SUCCESS );
			break;
		case 'w':
			sim->azi =- 2;
			break;
		case 's':
			sim->azi += 2;
			break;
		case 'd':
			sim->alt -= 2;
			break;
		case 'a':
			sim->alt += 2;
			break;
		case 'e':
			sim->zen -= 2;
			break;
		case 'q':
			sim->zen += 2;
			break;
		case 'z':
			sim->zoom -= 0.1;
			break;
		case 'c':
			sim->zoom += 0.1;
			break;
		case 't':
			if( sim->tracing ) {
				sim->tracing = 0;
				for( int i = 0; i < sim->nBodies; i ++ ) {
					for( int j = 0; j < NUMTRACE; j ++ ) {
						for( int k = 0; k < 3; k ++ ) {
							sim->bodies[ i ].trace[ j ][ k ] = 0;
						}
					}
				}
			}
			else if( !sim->tracing ) {
				sim->tracing = 1;
				for( int i = 0; i < sim->nBodies; i ++ ) {
					for( int j = 0; j < NUMTRACE; j ++ ) {
						for( int k = 0; k < 3; k ++ ) {
							sim->bodies[ i ].trace[ j ][ k ] = sim->bodies[ i ].pos[ k ];
						}
					}
				}
			}
			break;
		case 'r':
			printf( TEXT_YELLOW "Restarting: " TEXT_WHITE );
			zero( sim->nBodies );
			break;
		case 'b':
			zero( sim->nBodies += 3 );
			break;
		case 'v':
			if( sim->nBodies > 2 ) {
				zero( sim->nBodies -= 3 );
			}
			break;
		case 'k':
			save( );
			break;
		case 'l':
			load( );
			break;
	}
	return;
}

#define GRAVITY -0.5
void simulate( void ) {
	obj_t *this, *that;
	double tmp[ 3 ], d[ 3 ], r, dv[ 3 ], f, fv[ 3 ];
	for( int i = 0; i < sim->nBodies; i ++ ) {
		this = &sim->bodies[ i ];
		for( int k = 0; k < 3; k ++ ) {
			this->netforce[ k ] = 0;
		}
		for( int j = 0; j < sim->nBodies; j ++ ) {
			if( i == j ) { continue; }
			that = &sim->bodies[ j ];
			for( int k = 0; k < 3; k ++ ) {
				d[ k ] = this->pos[ k ] - that->pos[ k ];
			}
			r = sqrt( ( d[ 0 ] * d[ 0 ] ) + ( d[ 1 ] * d[ 1 ] ) + ( d[ 2 ] * d[ 2 ] ) );
			for( int k = 0; k < 3; k ++ ) {
				dv[ k ] = d[ k ] / r;
			}
			f = 0.00002 * ( GRAVITY * this->mass * that->mass ) / ( 1 + ( r * r ) );
			for( int k = 0; k < 3; k ++ ) {
				fv[ k ] = dv[ k ] * f;
				this->netforce[ k ] += fv[ k ];
			}
            if( r <= ( this->rad + that->rad ) ) {
                this->collided = 1;
                that->collided = 1;
                //compute final velocity from collisions
                for( int k = 0; k < 3; k ++ ) {
                    tmp[ k ] = this->vel[ k ];
                    this->vel[ k ] = 0.98 * ( ( ( this->mass - that->mass ) / ( this->mass + that->mass ) ) * this->vel[ k ] ) + ( ( ( 2 * that->mass ) / ( this->mass + that->mass ) ) * that->vel[ k ] );
                    this->vel[ k ] = 0.98 * ( ( ( 2 * this->mass ) / ( this->mass + that->mass ) ) * tmp[ k ] ) + ( ( ( that->mass - that->mass ) / ( this->mass + that->mass ) ) * this->vel[ k ] );
                }
            }
		}
		for( int k = 0; k < 3; k ++ ) {
			this->acc[ k ] = this->netforce[ k ] / this->mass;
            this->vel[ k ] = this->vel[ k ] + ( this->acc[ k ] * TIMESTEP );
            this->pos[ k ] = this->pos[ k ] + ( this->vel[ k ] * TIMESTEP );
		}
		if( sim->tracing ) {
	        for( int j = ( NUMTRACE - 1 ); j != 0; j -- ) {
	        	for( int k = 0; k < 3; k ++ ) {
	                this->trace[ j ][ k ] = this->trace[ j - 1 ][ k ];
	            }
	        }
	        //add current position to trace element 0
	        for( int k = 0; k < 3; k ++ ) {
	            this->trace[ 0 ][ k ] = this->pos[ k ];
	        }
	    }	
	}
}

void draw( void ) {
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glLoadIdentity( );
	glMatrixMode( GL_PROJECTION );
	glRotatef( sim->azi, 1.0, 0.0, 0.0 );
	glRotatef( sim->alt, 0.0, 1.0, 0.0 );
	glRotatef( sim->zen, 0.0, 0.0, 1.0 );
	glScalef( sim->zoom, sim->zoom, sim->zoom);
    for( int i = 0; i < sim->nBodies; i ++ ) {
    	obj_t *this = &sim->bodies[ i ];
        //draw traces
        if( sim->tracing ) {
            for( int j = 1; j < NUMTRACE; j ++ ) {
                //drawing lines with transparency is very computationally expensive, creates a noticable drop in frames per second when enabled
                glBegin( GL_LINES );
                    glColor4f( this->col[ 0 ], this->col[ 1 ], this->col[ 2 ], sim->alpha[ NUMTRACE - j ] );
                    glVertex3f( this->trace[ j - 1 ][ 0 ], this->trace[ j - 1 ][ 1 ], this->trace[ j - 1 ][ 2 ] );
                    glVertex3f( this->trace[ j ][ 0 ], this->trace[ j ][ 1 ], this->trace[ j ][ 2 ] );
                glEnd( );
            }
        }
        //draw individual bodies
        glPushMatrix( );
        glTranslatef( this->pos[ 0 ], this->pos[ 1 ], this->pos[ 2 ] );
        glColor3f( this->col[ 0 ], this->col[ 1 ], this->col[ 2 ] );
        glutSolidSphere( this->rad, 20, 20 );
        glPopMatrix( );
    }
    glutSwapBuffers( );
    glFlush( );
    return;
}

void save( void ) {
	FILE *output = fopen( "state", "w" );
	fprintf( output, "num = %i, time = %i\n", sim->nBodies, sim->time );
	for( int i = 0; i < sim->nBodies; i ++ ) {
		obj_t *this = &sim->bodies[ i ];
		fprintf( output, "pos = [ %lf %lf %lf ]\n", this->pos[ 0 ], this->pos[ 1 ], this->pos[ 2 ] );
		fprintf( output, "vel = [ %lf %lf %lf ]\n", this->vel[ 0 ], this->vel[ 1 ], this->vel[ 2 ] );
		fprintf( output, "acc = [ %lf %lf %lf ]\n", this->acc[ 0 ], this->acc[ 1 ], this->acc[ 2 ] );
		fprintf( output, "net = [ %lf %lf %lf ]\n", this->netforce[ 0 ], this->netforce[ 1 ], this->netforce[ 2 ] );
		fprintf( output, "rad = %lf, mass = %lf\n", this->rad, this->mass );
	}
	fclose( output );
	printf( TEXT_GREEN "Saved simulation to file.\n" TEXT_WHITE );
}

void load( void ) {
	FILE *input = fopen( "state", "r" );
	if( input == NULL ) { 
		printf( TEXT_RED "File not found.\n" TEXT_WHITE );
		return; }
	fscanf( input, "num = %i, time = %i\n", &sim->nBodies, &sim->time );
		for( int i = 0; i < sim->nBodies; i ++ ) {
		obj_t *this = &sim->bodies[ i ];
		fscanf( input, "pos = [ %lf %lf %lf ]\n", &this->pos[ 0 ], &this->pos[ 1 ], &this->pos[ 2 ] );
		fscanf( input, "vel = [ %lf %lf %lf ]\n", &this->vel[ 0 ], &this->vel[ 1 ], &this->vel[ 2 ] );
		fscanf( input, "acc = [ %lf %lf %lf ]\n", &this->acc[ 0 ], &this->acc[ 1 ], &this->acc[ 2 ] );
		fscanf( input, "net = [ %lf %lf %lf ]\n", &this->netforce[ 0 ], &this->netforce[ 1 ], &this->netforce[ 2 ] );
		fscanf( input, "rad = %lf, mass = %lf\n", &this->rad, &this->mass );
	}
	fclose( input );	
	printf( TEXT_GREEN "Loaded simulation from file.\n" TEXT_WHITE );
}

void timer( ) {
    draw( );
	simulate( );
	glutTimerFunc( TIMESTEP, timer, 0 );
	sim->time += TIMESTEP;
}

int main( int argc, char **argv ) {
	srand( time( NULL ) );
	zero( 2 );

	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowSize( WIDTH, HEIGHT );
	glutCreateWindow( "Gravity" );
	glEnable( GL_DEPTH_TEST );

	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glEnable( GL_BLEND );
	glClearColor( 0.0, 0.0, 0.0, 0.0 );

	timer( );
	glutKeyboardFunc( keys );
	glutDisplayFunc( draw );
	glutMainLoop( );

    return 0;
}
