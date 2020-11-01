//a basic, 3d gravity simulator written in C

/*  TODO
    -dynamic memory allocation when reading save file or arbritary length
    -implementing more accurate euclidian collision method
    -multithreading for read/write [save/load] functions
*/

//libaries
#include<stdio.h>
#include<GL/glut.h> //openGL
#include<stdlib.h>
#include<math.h> //math library for sqrt
#include<time.h>

//macros
#define WIDTH	 	750
#define HEIGHT 		750
#define NUMTRACE 	50
#define TIMESTEP 	20
#define GRAVITY     -0.1
#define TEXT_RED	"\033[31m"
#define TEXT_YELLOW "\033[33m"
#define TEXT_GREEN	"\033[32m"
#define TEXT_WHITE 	"\033[37m"

int i,j,k; //temporary counter variables

//planet structure
typedef struct
{
	double pos[3];
	double vel[3];
	double acc[3];
	double rad;
	double mass;
	double netforce[3];
	double trace[NUMTRACE][3];
	float col[3];
	_Bool collided;
	_Bool tracing;
} obj_t; //alias type

//universe structure
typedef struct
{
	obj_t *bodies;
	int nBodies;
	int time;
	float zoom,alpha[1000];
	double azi,alt,zen;
	_Bool tracing;
} uni_t; 

uni_t *sim; //alias type

//program function prototypes
double randDouble(double min,double max);
void init(int nBodies);
void save(void);
void load(void);
void simulate(void);
void draw(void);
void keys(unsigned char key,int x,int y);
void timer();

//generate random double in a given range
double randDouble(double min,double max)
{
    double random=((double)rand())/(double)RAND_MAX;
    double range=max-min;
    return(random*range)+min;
}

//zeroes out all bodies in simulation
void init(int nBodies) 
{
	sim=(uni_t*)malloc(sizeof(uni_t));
	sim->nBodies=nBodies;
	sim->bodies=(obj_t*)malloc(sizeof(obj_t)*sim->nBodies);
	sim->time=0;
	sim->zoom=1.0;
	for(i=0;i<NUMTRACE;i++){
		sim->alpha[i]=(float)(i+1)/NUMTRACE;
	}
	sim->azi=0.0;
	sim->alt=0.0;
	sim->zen=0.0;
	sim->tracing=0;
	obj_t *this=sim->bodies;
	for(i=0;i<sim->nBodies;i++){
		this->rad=randDouble(0.001,0.01);
		this->mass=4.19*this->rad; //approximation of volume/density/mass
		for(j=0;j<3;j++){
			this->pos[j]=randDouble(0.0,2.0)-1.0; //randomise position
			this->col[j]=randDouble(0.0,1.0); //randomise color
		}
		for(j=0;j<NUMTRACE;j++){
			for(k=0;k<3;k++){
				this->trace[j][k]=this->pos[k];
			}
		}
		this++; //increment body pointer
	}
	printf("Initialised simulation with %i bodies.\n",nBodies);
	return;
}

//get and handle keypress through openGL
void keys(unsigned char key,int x,int y) 
{
	int win;
	switch(key){
		case 27: //ESCAPE key
			win=glutGetWindow();
			glutDestroyWindow(win);
			exit(EXIT_SUCCESS);
			break;
		case 'w':
			sim->azi=-2;
			break;
		case 's':
			sim->azi+=2;
			break;
		case 'd':
			sim->alt-=2;
			break;
		case 'a':
			sim->alt+=2;
			break;
		case 'e':
			sim->zen-=2;
			break;
		case 'q':
			sim->zen+=2;
			break;
		case 'z':
			sim->zoom-=0.1;
			break;
		case 'c':
			sim->zoom+=0.1;
			break;
		case 't':
		    //turn off tracing by deleting all trace values
			if(sim->tracing){
				sim->tracing=0;
				for(i=0;i<sim->nBodies;i++){
					for(j=0;j<NUMTRACE;j++){
						for(k=0;k<3;k++){
							sim->bodies[i].trace[j][k]=0;
						}
					}
				}
			}
			//turn on tracing and create initial trace values
			else if(!sim->tracing){
				sim->tracing=1;
				for(i=0;i<sim->nBodies;i++){
					for(j=0;j<NUMTRACE;j++){
						for(k=0;k<3;k++){
							sim->bodies[i].trace[j][k]=sim->bodies[i].pos[k];
						}
					}
				}
			}
			break;
		case 'r':
			printf(TEXT_YELLOW"Restarting: "TEXT_WHITE);
			init(sim->nBodies);
			break;
		case 'b':
			init(sim->nBodies+=3);
			break;
		case 'v':
			if(sim->nBodies>2){ 
				init(sim->nBodies-=3);
			}
			break;
		case 'k':
			save();
			break;
		case 'l':
			load();
			break;
	}
	return;
}

//simulate between objects and gravity
void simulate(void)
{
    obj_t *this, *that; //*this and *that are structure pointers that are iterated for every body in simulation
    double tmp[3],d[3],r,dv[3],f,fv[3];
    for(i=0;i<sim->nBodies;i++){
        this=&sim->bodies[i]; //set this pointer to body[i]
        for(k=0;k<3;k++){
            this->netforce[k]=0; //set netforce of body[i] to zero
        }
        for(j=0;j<sim->nBodies;j++){ //iterate through all other bodies
            if(i==j){continue;} //dont compute on self
            that=&sim->bodies[j]; //set that pointer to body[j]
            for(k=0;k<3;k++){
                d[k]=this->pos[k]-that->pos[k]; //distance between body[i] and body[j]
            }
            r=sqrt((d[0]*d[0])+(d[1]*d[1])+(d[2]*d[2])); //magnitude of d[k]
            for(k=0;k<3;k++){
                dv[k]=d[k]/r; //normalised distance
            }
            f=0.00002*(GRAVITY*this->mass*that->mass)/(1+(r*r)); //compute force due to gravity
            for(k=0;k<3;k++){
                fv[k]=dv[k]*f; //normalised force
                this->netforce[k]+=fv[k];
            }
            if(r<=(this->rad+that->rad)){ //check to see if bodies have collided
                this->collided=1;
                that->collided=1;
                //compute reflected velocity from collisions
                for(k=0;k<3;k++){
                    tmp[k]=this->vel[k];
                    //a quick and dirty method of computing euclidian reflections in 3d
                    this->vel[k]=0.98*(((this->mass-that->mass)/(this->mass+that->mass))*this->vel[k])+(((2*that->mass)/(this->mass+that->mass))*that->vel[k]);
                }
            }
		}
		for(k=0;k<3;k++){
			this->acc[k]=this->netforce[k]/this->mass; //update acceleration from netforce
            this->vel[k]=this->vel[k]+(this->acc[k]*TIMESTEP); //update velocity from acceleration
            this->pos[k]=this->pos[k]+(this->vel[k]*TIMESTEP); //update position from velocity
		}
		if(sim->tracing) { //check if tracing
	        for(j=(NUMTRACE-1);j!=0;j--){ //loop backwards through trace positions
	        	for(k=0;k<3;k++){
	                this->trace[j][k]=this->trace[j-1][k]; //draw traces as required
	            }
	        }
	        //add current position to trace element 0
	        for(k=0;k<3;k++){
	            this->trace[0][k]=this->pos[k];
	        }
	    }	
	}  
}

//draw game data to screen - openGL primitives
void draw(void)
{
    //setup view in space
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glRotatef(sim->azi,1.0,0.0,0.0);
	glRotatef(sim->alt,0.0,1.0,0.0);
	glRotatef(sim->zen,0.0,0.0,1.0);
	glScalef(sim->zoom,sim->zoom,sim->zoom);
	//begin looping through bodies
    for(i=0;i<sim->nBodies;i++){
    	obj_t *this=&sim->bodies[i];
        //draw traces
        if(sim->tracing){
            for(j=1;j<NUMTRACE;j++){
                //drawing lines with transparency is computationally expensive, creates a noticable drop in frames per second when enabled
                glBegin(GL_LINES);
                    glColor4f(this->col[0],this->col[1],this->col[2],sim->alpha[NUMTRACE-j]); //this will enable alpha tracing, which as mentioned is computationally expensive
                    glVertex3f(this->trace[j-1][0],this->trace[j-1][1],this->trace[j-1][2]);
                    glVertex3f(this->trace[j][0],this->trace[j][1],this->trace[j][2]);
                glEnd();
            }
        }
        //draw individual bodies
        glPushMatrix();
        glTranslatef(this->pos[0],this->pos[1],this->pos[2]);
        glColor3f(this->col[0],this->col[1],this->col[2]);
        glutSolidSphere(this->rad,20,20);
        glPopMatrix();
    }
    //refresh display and return
    glutSwapBuffers();
    glFlush();
    return;
}

//write game data to text file
void save(void) 
{
	FILE *output=fopen("state","w"); //opens or creates new file
	fprintf(output,"num = %i, time = %i\n",sim->nBodies,sim->time); //print header
    //loop through all bodies and serialise data
	for(i=0;i<sim->nBodies;i++){
		obj_t *this=&sim->bodies[i];
		fprintf(output,"pos = [ %lf %lf %lf ]\n",this->pos[0],this->pos[1],this->pos[2]);
		fprintf(output,"vel = [ %lf %lf %lf ]\n",this->vel[0],this->vel[1],this->vel[2]);
		fprintf(output,"acc = [ %lf %lf %lf ]\n",this->acc[0],this->acc[1],this->acc[2]);
		fprintf(output,"net = [ %lf %lf %lf ]\n",this->netforce[0],this->netforce[1],this->netforce[2]);
		fprintf(output,"rad = %lf, mass = %lf\n",this->rad,this->mass);
	}
	//close file
	fclose(output);
	printf(TEXT_GREEN"Saved simulation to file.\n"TEXT_WHITE);
}

//load in game data from text file
void load(void)
{
	FILE *input=fopen("state","r");
	//check if file exists, if not return
	if(input==NULL){ 
		printf(TEXT_RED"File not found.\n"TEXT_WHITE);
		return; 
	}
	//stream file data into structure
	fscanf(input,"num = %i, time = %i\n",&sim->nBodies,&sim->time);
		for(i=0;i<sim->nBodies;i++){
		obj_t *this=&sim->bodies[i];
		fscanf(input,"pos = [ %lf %lf %lf ]\n",&this->pos[0],&this->pos[1],&this->pos[2]);
		fscanf(input,"vel = [ %lf %lf %lf ]\n",&this->vel[0],&this->vel[1],&this->vel[2]);
		fscanf(input,"acc = [ %lf %lf %lf ]\n",&this->acc[0],&this->acc[1],&this->acc[2]);
		fscanf(input,"net = [ %lf %lf %lf ]\n",&this->netforce[0],&this->netforce[1],&this->netforce[2]);
		fscanf(input,"rad = %lf, mass = %lf\n",&this->rad,&this->mass);
	}
	fclose(input);	
	printf(TEXT_GREEN"Loaded simulation from file.\n"TEXT_WHITE);
}

//timer function - calls critical program code
void timer() //this also takes no variables, not even 'void' 
{
    draw();
	simulate();
	glutTimerFunc(TIMESTEP,timer,0);
	sim->time+=TIMESTEP;
}

//main function - initialise and then handoff to openGLs
int main(int argc,char **argv)
{
	srand(time(NULL));
	init(2);

	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutInitWindowSize(WIDTH,HEIGHT);
	glutCreateWindow("Gravity");
	glEnable(GL_DEPTH_TEST);

	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glClearColor(0.0,0.0,0.0,0.0);

	timer();
	glutKeyboardFunc(keys);
	glutDisplayFunc(draw);
	glutMainLoop();

    return 0;
}
