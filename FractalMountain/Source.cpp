#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <random>
#include <math.h>
#include <GL/glut.h>
#include <time.h>
#define PI 3.14159265

bool drawFogFlag = false;



int lowQualityDepth = 7;
int mediumQualityDepth = 8;
int highQualityDepth = 9;
int ultraHighQualityDepth = 10;

const char * daySkyboxFilename = "./day.bmp";
const char * nightSkyboxFilename = "./night.bmp";
const char * eveningSkyboxFilename = "./evening.bmp";
const char * interstellarSkyboxFilename = "./interstellar.bmp";

const char * seaFilename = "./sea_seamless.bmp";
float seaHeightNormal = 6 / 16.0;

const char * elevationFilename = "./elevation.bmp";




double lowRoughnessD = 2.1f;
double lowRoughnessR = 0.2f;
double lowRoughnessSmoothihgPasses = 3;

double mediumRoughnessD = 2.2f;
double mediumRoughnessR = 0.3f;
double mediumRoughnessSmoothihgPasses = 2;

double highRoughnessD = 2.3f;
double highRoughnessR = 0.4f;
double highRoughnessSmoothihgPasses = 1;



GLuint skyboxTexture;
bool drawSkyboxFlag = true;


GLuint elevationMap;
bool useElevationMap = true;


GLuint seaTexture;
GLuint drawSeaFlag = true;



typedef GLfloat point3[3];
typedef GLfloat point2[2];


point3** terrain = NULL;
float terrainSize = 50;

int terrainDepth = mediumQualityDepth;
double D = mediumRoughnessD;
double R = mediumRoughnessR;
double smoothingPasses = mediumRoughnessSmoothihgPasses;

bool drawWireframeFlag = false;
bool drawHeightmapFlag = false;


float cameraRadius = 110;
point3 cameraPosition = { 0, 45, cameraRadius };
point3 cameraCenter = { 0, 45, 0 };
float cameraAngle = 0;
float cameraRotationPerSec = 60.0;
GLfloat fieldOfView = 60.0f;

float lightRadius = 60;
point3 lightPosition = { lightRadius, 0, 0 };
point3 lightCenter{ 0, 0, 0 };
float lightAngle = 0;
float lightRotationPerSec = 30;



long previousTime;
long currentTime;

bool keyLeft = false;
bool keyRight = false;

bool leftButton = false;
bool rightButton = false;


long displayPrev;

long displayCurrent;


point3 v[] = { { 0.0, 0.0, 1.0 }, { 0.0, 0.942809, -0.33333 }, { -0.816497, -0.471405, -0.333333 }, { 0.816497, -0.471405, -0.333333 } };



GLuint loadBMP_custom(const char * imagepath){
	// Data read from the header of the BMP file
	unsigned char header[54]; // Each BMP file begins by a 54-bytes header
	unsigned int dataPos;     // Position in the file where the actual data begins
	unsigned int width, height;
	unsigned int imageSize;   // = width*height*3
	// Actual RGB data
	unsigned char * data;

	FILE * file = fopen(imagepath, "rb");
	if (!file){
		printf("Image could not be opened\n");
		return 0;
	}

	// If not 54 bytes read : problem
	if (fread(header, 1, 54, file) != 54){
		printf("Not a correct BMP file\n");
		return false;
	}

	if (header[0] != 'B' || header[1] != 'M'){
		printf("Not a correct BMP file\n");
		return 0;
	}


	// Read ints from the byte array
	dataPos = *(int*)&(header[0x0A]);
	imageSize = *(int*)&(header[0x22]);
	width = *(int*)&(header[0x12]);
	height = *(int*)&(header[0x16]);



	// Some BMP files are misformatted, guess missing information
	if (imageSize == 0)    imageSize = width*height * 3; // 3 : one byte for each Red, Green and Blue component
	if (dataPos == 0)      dataPos = 54; // The BMP header is done that way


	// Create a buffer
	data = new unsigned char[imageSize];

	// Read the actual data from the file into the buffer
	fread(data, 1, imageSize, file);

	//Everything is in memory now, the file can be closed
	fclose(file);

	if (false){
		for (int i = 0; i < imageSize; i += 3){
			int tmp = data[i];
			data[i] = data[i + 2];
			data[i + 2] = tmp;
			std::cout << i / 3 << " -- R:" << (int)data[i] << " G:" << (int)data[i + 1] << " B:" << (int)data[i + 2] << std::endl << std::endl;
		}
	}

	// Create one OpenGL texture
	GLuint textureID;
	glGenTextures(1, &textureID);

	// "Bind" the newly created texture : all future texture functions will modify this texture
	glBindTexture(GL_TEXTURE_2D, textureID);

	// Give the image to OpenGL
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, data);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	return textureID;
}




void drawSkybox()
{
	glBindTexture(GL_TEXTURE_2D, skyboxTexture);


	// Store the current matrix
	glPushMatrix();

	// Reset and transform the matrix.
	glLoadIdentity();
	gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2],
		0, 0, 0,
		0, 1, 0);

	// Enable/Disable features
	glPushAttrib(GL_ENABLE_BIT);

	glDisable(GL_FOG);
	glDisable(GL_LIGHTING);

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);

	// Just in case we set all vertices to white.

	float skyboxDistance = 600;

	// Render the front quad
	glBegin(GL_QUADS);
	glTexCoord2f(1 / 4.0, 1 / 3.0); glVertex3f(skyboxDistance, -skyboxDistance, -skyboxDistance);
	glTexCoord2f(2 / 4.0, 1 / 3.0); glVertex3f(-skyboxDistance, -skyboxDistance, -skyboxDistance);
	glTexCoord2f(2 / 4.0, 2 / 3.0); glVertex3f(-skyboxDistance, skyboxDistance, -skyboxDistance);
	glTexCoord2f(1 / 4.0, 2 / 3.0); glVertex3f(skyboxDistance, skyboxDistance, -skyboxDistance);
	glEnd();

	// Render the left quad
	glBegin(GL_QUADS);
	glTexCoord2f(0, 1 / 3.0); glVertex3f(skyboxDistance, -skyboxDistance, skyboxDistance);
	glTexCoord2f(1 / 4.0, 1 / 3.0); glVertex3f(skyboxDistance, -skyboxDistance, -skyboxDistance);
	glTexCoord2f(1 / 4.0, 2 / 3.0); glVertex3f(skyboxDistance, skyboxDistance, -skyboxDistance);
	glTexCoord2f(0, 2 / 3.0); glVertex3f(skyboxDistance, skyboxDistance, skyboxDistance);
	glEnd();

	// Render the back quad
	glBegin(GL_QUADS);
	glTexCoord2f(3 / 4.0, 1 / 3.0); glVertex3f(-skyboxDistance, -skyboxDistance, skyboxDistance);
	glTexCoord2f(4 / 4.0, 1 / 3.0); glVertex3f(skyboxDistance, -skyboxDistance, skyboxDistance);
	glTexCoord2f(4 / 4.0, 2 / 3.0); glVertex3f(skyboxDistance, skyboxDistance, skyboxDistance);
	glTexCoord2f(3 / 4.0, 2 / 3.0); glVertex3f(-skyboxDistance, skyboxDistance, skyboxDistance);

	glEnd();

	// Render the right quad
	glBegin(GL_QUADS);
	glTexCoord2f(2 / 4.0, 1 / 3.0); glVertex3f(-skyboxDistance, -skyboxDistance, -skyboxDistance);
	glTexCoord2f(3 / 4.0, 1 / 3.0); glVertex3f(-skyboxDistance, -skyboxDistance, skyboxDistance);
	glTexCoord2f(3 / 4.0, 2 / 3.0); glVertex3f(-skyboxDistance, skyboxDistance, skyboxDistance);
	glTexCoord2f(2 / 4.0, 2 / 3.0); glVertex3f(-skyboxDistance, skyboxDistance, -skyboxDistance);
	glEnd();

	// Render the top quad
	glBegin(GL_QUADS);
	glTexCoord2f(2 / 4.0, 2 / 3.0); glVertex3f(-skyboxDistance, skyboxDistance, -skyboxDistance);
	glTexCoord2f(2 / 4.0, 1); glVertex3f(-skyboxDistance, skyboxDistance, skyboxDistance);
	glTexCoord2f(1 / 4.0, 1); glVertex3f(skyboxDistance, skyboxDistance, skyboxDistance);
	glTexCoord2f(1 / 4.0, 2 / 3.0); glVertex3f(skyboxDistance, skyboxDistance, -skyboxDistance);
	glEnd();


	// Render the bottom quad
	glBegin(GL_QUADS);
	glTexCoord2f(2 / 4.0, 1 / 3.0); glVertex3f(-skyboxDistance, -skyboxDistance, -skyboxDistance);
	glTexCoord2f(1 / 4.0, 1 / 3.0); glVertex3f(skyboxDistance, -skyboxDistance, -skyboxDistance);
	glTexCoord2f(1 / 4.0, 0); glVertex3f(skyboxDistance, -skyboxDistance, skyboxDistance);
	glTexCoord2f(2 / 4.0, 0);  glVertex3f(-skyboxDistance, -skyboxDistance, skyboxDistance);
	glEnd();

	// Restore enable bits and matrix
	glPopAttrib();
	glPopMatrix();


}


void calculateTerrain(point3 a, point3 b, point3 c, point3 d)
{
	int terrainDims = pow(2, terrainDepth) + 1;


	terrain = new point3*[terrainDims];
	for (int i = 0; i < terrainDims; i++){
		terrain[i] = new point3[terrainDims];
	}

	//Set the terrains x and z values for each point
	for (int i = 0; i < terrainDims; i++){
		for (int j = 0; j < terrainDims; j++){
			terrain[i][j][0] = a[0] + j * (b[0] - a[0]) / (terrainDims - 1);
			terrain[i][j][2] = a[2] + i * (d[2] - a[2]) / (terrainDims - 1);
		}
	}
	//Set the y value for the initial points
	terrain[0][0][1] = a[1];
	terrain[0][terrainDims - 1][1] = b[1];
	terrain[terrainDims - 1][terrainDims - 1][1] = c[1];
	terrain[terrainDims - 1][0][1] = d[1];


	std::random_device rd;
	std::mt19937 gen(rd());



	for (int k = terrainDepth; k > 0; k--){

		double H = 3 - D;

		double squareSize = (b[0] - a[0]) / pow(2, terrainDepth - k);
		if (squareSize < 0) squareSize = -squareSize;

		float tmp = pow(squareSize, 2 * (3 - D));
		//float tmp = (1 - pow(2, 2 * H - 2)); 
		//float tmp = pow(pow(2, terrainDepth - k), 2 * H); //

		tmp = R*sqrt(tmp);
		//tmp = 0.2*squareSize;


		std::normal_distribution<double> normalDistribution(0, tmp);


		int squareSizeInGrid = pow(2, k);

		for (int i = 0; i < terrainDims - 1; i += squareSizeInGrid){
			for (int j = 0; j < terrainDims - 1; j += squareSizeInGrid){

				double error;

				//if first row, calculate top
				if (i == 0){
					error = normalDistribution(gen);//+(1.0 *rand() / RAND_MAX) * 2 * tmp - tmp;
					terrain[0][j + squareSizeInGrid / 2][1] = (terrain[0][j][1] + terrain[0][j + squareSizeInGrid][1]) / 2 + error;
				}

				//if first column, calculate left
				if (j == 0){
					error = normalDistribution(gen);
					terrain[i + squareSizeInGrid / 2][0][1] = (terrain[i][0][1] + terrain[i + squareSizeInGrid][0][1]) / 2 + error;
				}

				//calculate right
				error = normalDistribution(gen);
				terrain[i + squareSizeInGrid / 2][j + squareSizeInGrid][1] = (terrain[i][j + squareSizeInGrid][1] + terrain[i + squareSizeInGrid][j + squareSizeInGrid][1]) / 2 + error;

				//calculate bottom
				error = normalDistribution(gen);
				terrain[i + squareSizeInGrid][j + squareSizeInGrid / 2][1] = (terrain[i + squareSizeInGrid][j][1] + terrain[i + squareSizeInGrid][j + squareSizeInGrid][1]) / 2 + error;

				//calculate center
				error = normalDistribution(gen);
				terrain[i + squareSizeInGrid / 2][j + squareSizeInGrid / 2][1] = (terrain[i][j][1] + terrain[i][j + squareSizeInGrid][1] + terrain[i + squareSizeInGrid][j + squareSizeInGrid][1] + terrain[i + squareSizeInGrid][j][1]) / 4 + error;
			}
		}
	}

	//Smoothing Filter
	float** smoothedHeight = new float*[terrainDims];
	for (int i = 0; i < terrainDims; i++){
		smoothedHeight[i] = new float[terrainDims];
	}

	float smoothness = 0.8;

	for (int k = 0; k < smoothingPasses; k++){

		for (int i = 0; i < terrainDims; i++){
			for (int j = 0; j < terrainDims; j++){
				float temp = 0;
				int noNeighboors = 0;
				//up
				if (i > 0){
					temp += terrain[i - 1][j][1];
					noNeighboors++;
				}

				//up - left
				if (i > 0 && j > 0){
					temp += terrain[i - 1][j - 1][1];
					noNeighboors++;
				}

				//up - left
				if (i > 0 && j < terrainDims - 1){
					temp += terrain[i - 1][j + 1][1];
					noNeighboors++;
				}



				//down
				if (i < terrainDims - 1){
					temp += terrain[i + 1][j][1];
					noNeighboors++;
				}



				//down - left
				if (i < terrainDims - 1 && j > 0){
					temp += terrain[i + 1][j - 1][1];
					noNeighboors++;
				}


				//down - right
				if (i < terrainDims - 1 && j < terrainDims - 1){
					temp += terrain[i + 1][j + 1][1];
					noNeighboors++;
				}


				//left
				if (j > 0){
					temp += terrain[i][j - 1][1];
					noNeighboors++;
				}


				//right
				if (j < terrainDims - 1){
					temp += terrain[i][j + 1][1];
					noNeighboors++;
				}


				temp /= noNeighboors;

				smoothedHeight[i][j] = temp + (1 - smoothness) * (terrain[i][j][1] - temp);
				//terrain[i][j][1] = temp + (1 - smoothness) * (terrain[i][j][1] - temp);
			}
		}


		for (int i = 0; i < terrainDims; i++){
			for (int j = 0; j < terrainDims; j++){
				terrain[i][j][1] = smoothedHeight[i][j];
			}
		}
	}


	for (int i = 0; i < terrainDims; i++){
		free(smoothedHeight[i]);
	}
	free(smoothedHeight);


}



void drawHeightmap(){

	glLoadIdentity();
	gluLookAt(0, 0, 0,
		0, 0, -1,
		0, 1, 0);




	glPushAttrib(GL_ENABLE_BIT);

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	//glDisable(GL_DEPTH_TEST);
	//glDisable(GL_NORMALIZE);
	//glDisable(GL_LINEAR_ATTENUATION);
	glDisable(GL_FOG);
	glEnable(GL_BLEND);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-50.0, 50.0, -50.0, 50.0);
	glMatrixMode(GL_MODELVIEW);

	int terrainDims = pow(2, terrainDepth) + 1;

	float max = terrain[0][0][1];
	float min = terrain[0][0][1];

	for (int i = 0; i < terrainDims - 1; i++){
		for (int j = 0; j < terrainDims - 1; j++){
			if (terrain[i][j][1] > max){
				max = terrain[i][j][1];
			}
			if (terrain[i][j][1] < min){
				min = terrain[i][j][1];
			}

		}
	}



	for (int i = 0; i < terrainDims - 1; i++){
		for (int j = 0; j < terrainDims - 1; j++){


			float width = 100 / pow(2, terrainDepth);

			point2 a;
			a[0] = -50 + j*width;
			a[1] = 50 - i*width;

			point2 b;
			b[0] = a[0] + width;
			b[1] = a[1];

			point2 c;
			c[0] = a[0] + width;
			c[1] = a[1] - width;

			point2 d;
			d[0] = a[0];
			d[1] = a[1] - width;


			float col = (terrain[i][j][1] - min) / (max - min);

			glBegin(GL_QUADS);
			glColor3f(col, col, col);
			glVertex2fv(a);
			glVertex2fv(b);
			glVertex2fv(c);
			glVertex2fv(d);

			glEnd();

		}
	}

	glPopAttrib();
}


void getNormal(point3 normal, point3 a, point3 b, point3 c){
	//AB
	point3 A = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
	//AC
	point3 B = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };

	normal[0] = A[1] * B[2] - A[2] * B[1];
	normal[1] = A[2] * B[0] - A[0] * B[2];
	normal[2] = A[0] * B[1] - A[1] * B[0];
}


void drawTerrain()
{
	float zeros[3] = { 0, 0, 0 };
	GLfloat whiteMaterial[] = { 1, 1, 1 };



	glMaterialfv(GL_FRONT, GL_SPECULAR, zeros);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, whiteMaterial);
	glMaterialfv(GL_FRONT, GL_AMBIENT, whiteMaterial);


	glPushAttrib(GL_ENABLE_BIT);
	glPushMatrix();
	if (!drawWireframeFlag){
		glEnable(GL_TEXTURE_2D);
		glEnable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
		//glCullFace(GL_BACK);
		//glEnable(GL_CULL_FACE);
		glBindTexture(GL_TEXTURE_2D, elevationMap);
	}
	else{

		glDisable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_NORMALIZE);
		glEnable(GL_LINEAR_ATTENUATION);
		glShadeModel(GL_SMOOTH);
	}



	int terrainDims = pow(2, terrainDepth) + 1;


	GLfloat max = terrain[0][0][1];
	GLfloat min = terrain[0][0][1];

	for (int i = 0; i < terrainDims - 1; i++){
		for (int j = 0; j < terrainDims - 1; j++){
			if (terrain[i][j][1] > max){
				max = terrain[i][j][1];
			}
			if (terrain[i][j][1] < min){
				min = terrain[i][j][1];
			}

		}
	}

	GLfloat height = max - min;


	for (int i = 0; i < terrainDims - 1; i++){
		for (int j = 0; j < terrainDims - 1; j++){

			point3 a;
			a[0] = terrain[i][j][0];
			a[1] = terrain[i][j][1];
			a[2] = terrain[i][j][2];

			point3 b;
			b[0] = terrain[i][j + 1][0];
			b[1] = terrain[i][j + 1][1];
			b[2] = terrain[i][j + 1][2];

			point3 c;
			c[0] = terrain[i + 1][j + 1][0];
			c[1] = terrain[i + 1][j + 1][1];
			c[2] = terrain[i + 1][j + 1][2];

			point3 d;
			d[0] = terrain[i + 1][j][0];
			d[1] = terrain[i + 1][j][1];
			d[2] = terrain[i + 1][j][2];


			//AB
			point3 A = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
			//AC
			point3 B = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };





			if (drawWireframeFlag){

				//Draw using quads
				point3 up;
				getNormal(up, a, b, c);
				glNormal3fv(up);
				glBegin(GL_LINE_STRIP);
				glVertex3fv(a);
				glVertex3fv(b);
				glVertex3fv(c);
				glVertex3fv(d);
				glEnd();

			}
			else{


				GLfloat aNorm = (a[1] - min) / height;
				if (aNorm < 0) aNorm = 0;
				if (aNorm > 1) aNorm = 1;
				GLfloat bNorm = (b[1] - min) / height;
				if (bNorm < 0) bNorm = 0;
				if (bNorm > 1) bNorm = 1;

				GLfloat cNorm = (c[1] - min) / height;
				if (cNorm < 0) cNorm = 0;
				if (cNorm > 1) cNorm = 1;

				GLfloat dNorm = (d[1] - min) / height;
				if (dNorm < 0) dNorm = 0;
				if (dNorm > 1) dNorm = 1;



				glBegin(GL_TRIANGLES);

				point3 up1 = { 0, 1, 0 };
				getNormal(up1, a, b, c);
				glNormal3fv(up1);

				glTexCoord2f(0.5, aNorm); glVertex3fv(a);
				glTexCoord2f(0.5, bNorm); glVertex3fv(b);
				glTexCoord2f(0.5, cNorm); glVertex3fv(c);
				glEnd();


				glBegin(GL_TRIANGLES);

				point3 up2 = { 0, 1, 0 };
				getNormal(up2, a, c, d);
				glNormal3fv(up2);

				glTexCoord2f(0.5, aNorm); glVertex3fv(a);
				glTexCoord2f(0.5, cNorm); glVertex3fv(c);
				glTexCoord2f(0.5, dNorm); glVertex3fv(d);

				glEnd();
			}

		}
	}
	glPopAttrib();
	glPopMatrix();
}


void drawSea(){


	float seaHeight;

	int terrainDims = pow(2, terrainDepth) + 1;


	GLfloat max = terrain[0][0][1];
	GLfloat min = terrain[0][0][1];

	for (int i = 0; i < terrainDims - 1; i++){
		for (int j = 0; j < terrainDims - 1; j++){
			if (terrain[i][j][1] > max){
				max = terrain[i][j][1];
			}
			if (terrain[i][j][1] < min){
				min = terrain[i][j][1];
			}

		}
	}

	float terrainHeight = max - min;

	seaHeight = min + seaHeightNormal*(terrainHeight);


	glPushAttrib(GL_ENABLE_BIT);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	glBindTexture(GL_TEXTURE_2D, seaTexture);

	glPushMatrix();
	glBegin(GL_QUADS);

	glTexCoord2f(0, 0);
	glVertex3f(-terrainSize, seaHeight, -terrainSize);
	glTexCoord2f(2, 0);
	glVertex3f(terrainSize, seaHeight, -terrainSize);
	glTexCoord2f(2, 2);
	glVertex3f(terrainSize, seaHeight, terrainSize);
	glTexCoord2f(0, 2);
	glVertex3f(-terrainSize, seaHeight, terrainSize);

	glEnd();

	glPopMatrix();
	glPopAttrib();

}



void normalize(point3 p)
{
	double sum = 0.0;
	sum += p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
	sum = sqrt(sum);
	if (sum > 0.0){
		p[0] /= sum;
		p[1] /= sum;
		p[2] /= sum;
	}
}




void triangle(point3 a, point3 b, point3 c)
{
	int mode = 0;

	glBegin(GL_LINE_LOOP);

	glVertex3fv(a);

	glVertex3fv(b);

	glVertex3fv(c);

	glEnd();
}



void divide_triangle(point3 a, point3 b, point3 c, int depth)
{
	point3 v1, v2, v3;
	int j;
	if (depth == 0){
		triangle(a, b, c); /* draw triangle at end of recursion */
		return;
	}


	for (j = 0; j < 3; j++) v1[j] = a[j] + b[j];
	normalize(v1);
	for (j = 0; j < 3; j++) v2[j] = a[j] + c[j];
	normalize(v2);
	for (j = 0; j < 3; j++) v3[j] = b[j] + c[j];
	normalize(v3);
	divide_triangle(a, v1, v2, depth - 1);
	divide_triangle(c, v2, v3, depth - 1);
	divide_triangle(b, v3, v1, depth - 1);
	divide_triangle(v1, v3, v2, depth - 1);
}



void tetrahedron(int m)
/* Apply triangle subdivision to faces of tetrahedron */
{
	divide_triangle(v[0], v[1], v[2], m);
	divide_triangle(v[3], v[2], v[1], m);
	divide_triangle(v[0], v[3], v[1], m);
	divide_triangle(v[0], v[2], v[3], m);
}






void drawSun(){

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_FOG);
	glEnable(GL_DEPTH_TEST);

	//Draw Sun
	GLfloat yellowMaterial[] = { 1, 1, 0 };
	float zeros[3] = { 0, 0, 0 };

	glMaterialfv(GL_FRONT, GL_SPECULAR, zeros);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, zeros);
	glMaterialfv(GL_FRONT, GL_AMBIENT, zeros);
	glMaterialfv(GL_FRONT, GL_EMISSION, yellowMaterial);

	glPushMatrix();
	glTranslatef(lightPosition[0], lightPosition[1], lightPosition[2]);
	glScalef(3, 3, 3);
	tetrahedron(6);
	glPopMatrix();



	glPopAttrib();
}








void moveCamera(float degrees){

	cameraAngle += degrees;

	//Vale tin gwnia sto [0,360)
	if (cameraAngle >= 360){
		cameraAngle -= 360;
	}
	if (cameraAngle < 0){
		cameraAngle += 360;
	};

	//cos = x/r, sin = z/r
	cameraPosition[0] = cameraCenter[0] + cameraRadius * cos(cameraAngle * PI / 180);
	cameraPosition[2] = cameraCenter[2] + cameraRadius * sin(cameraAngle * PI / 180);
}


void moveLight(float degrees){
	lightAngle += degrees;

	//Vale tin gwnia sto [0,180) kai otan ftasei stis 0 moires ksanapaei stis 180
	if (lightAngle <= 0){
		lightAngle = 180;
	}

	//cos = x/r, sin = z/r
	lightPosition[0] = lightCenter[0] + lightRadius * cos(lightAngle * PI / 180);
	lightPosition[1] = lightCenter[2] + lightRadius * sin(lightAngle * PI / 180);
}




void idle(){
	previousTime = currentTime;
	currentTime = clock();

	float timeInterval = (double)(currentTime - previousTime) / CLOCKS_PER_SEC;

	if (keyLeft || leftButton){
		moveCamera(timeInterval * cameraRotationPerSec);
	}

	if (keyRight || rightButton){
		moveCamera(-timeInterval * cameraRotationPerSec);
	}

	moveLight(-timeInterval * lightRotationPerSec);

	glutPostRedisplay();
}





void mouseFunction(int button, int state, int x, int y){

	//Ama pati8ike to koumpi
	if (state == GLUT_DOWN){
		if (button == GLUT_LEFT_BUTTON){
			leftButton = true;
		}
		if (button == GLUT_RIGHT_BUTTON){
			rightButton = true;
		}
	}
	//An to afise
	else if (state == GLUT_UP){
		if (button == GLUT_LEFT_BUTTON){
			leftButton = false;
		}
		if (button == GLUT_RIGHT_BUTTON){
			rightButton = false;
		}
	}
}



void keyPressed(unsigned char key, int x, int y) {
	if (key == 'a'){
		keyLeft = true;
	}
	if (key == 'd'){
		keyRight = true;
	}
}

void keyReleased(unsigned char key, int x, int y) {
	if (key == 'a'){
		keyLeft = false;
	}
	if (key == 'd'){
		keyRight = false;
	}
}

void specialKeyPressed(int key, int x, int y){
	if (key == GLUT_KEY_LEFT){
		keyLeft = true;
	}
	if (key == GLUT_KEY_RIGHT){
		keyRight = true;
	}
}

void specialKeyReleased(int key, int x, int y){
	if (key == GLUT_KEY_LEFT){
		keyLeft = false;
	}
	if (key == GLUT_KEY_RIGHT){
		keyRight = false;
	}
}





void RightClickMenu(int id){
	int numberOfPoints;

	//mode = id;
	if (id >= 1 && id <= 5){
		drawSkyboxFlag = true;
		drawHeightmapFlag = false;
	}

	if (id >= 6 && id <= 9){
		int terrainDims = pow(2, terrainDepth) + 1;
		if (terrain != NULL){
			for (int i = 0; i < terrainDims; i++){
				free(terrain[i]);
			}
			free(terrain);
		}
	}
	switch (id){
	case 1:
		skyboxTexture = loadBMP_custom(daySkyboxFilename);
		break;
	case 2:
		skyboxTexture = loadBMP_custom(nightSkyboxFilename);
		break;
	case 3:
		skyboxTexture = loadBMP_custom(eveningSkyboxFilename);
		break;
	case 4:
		skyboxTexture = loadBMP_custom(interstellarSkyboxFilename);
		break;
	case 5:
		drawSkyboxFlag = false;
		break;
	case 6:

		terrainDepth = lowQualityDepth;
		break;
	case 7:


		terrainDepth = mediumQualityDepth;
		break;
	case 8:


		terrainDepth = highQualityDepth;
		break;
	case 9:


		terrainDepth = ultraHighQualityDepth;
		break;
	case 10:
		D = lowRoughnessD;
		R = lowRoughnessR;
		smoothingPasses = lowRoughnessSmoothihgPasses;
		break;
	case 11:
		D = mediumRoughnessD;
		R = mediumRoughnessR;
		smoothingPasses = mediumRoughnessSmoothihgPasses;
		break;
	case 12:

		D = highRoughnessD;
		R = highRoughnessR;
		smoothingPasses = highRoughnessSmoothihgPasses;
		break;
	case 13:
		drawSeaFlag = !drawSeaFlag;
		drawHeightmapFlag = false;
		break;
	case 14:
		drawWireframeFlag = !drawWireframeFlag;
		drawHeightmapFlag = false;
		break;
	case 15:
		drawFogFlag = !drawFogFlag;
		drawHeightmapFlag = false;
		break;
	case 16:
		drawHeightmapFlag = !drawHeightmapFlag;
		break;
	default:
		exit(0);
		break;
	}

	if (!drawHeightmapFlag){

		int width = glutGet(GLUT_WINDOW_WIDTH);
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		glViewport(0, 0, (GLsizei)width, (GLsizei)height);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(fieldOfView, (GLfloat)width / (GLfloat)height, 1, 1500.0);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}



	if (id >= 6 && id <= 12){
		point3  a = { -terrainSize, 0, terrainSize };
		point3 b = { terrainSize, 0, terrainSize };
		point3 c = { terrainSize, 0, -terrainSize };
		point3 d = { -terrainSize, 0, -terrainSize };
		calculateTerrain(a, b, c, d);
	}


	glutPostRedisplay();
}








void display()
{
	displayPrev = displayCurrent;
	displayCurrent = clock();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	if (drawHeightmapFlag){

		drawHeightmap();
		glutSwapBuffers();
		return;
	}


	if (drawFogFlag){
		glEnable(GL_FOG);
		glFogf(GL_FOG_MODE, GL_EXP);
		glFogf(GL_FOG_DENSITY, 0.008);
		GLfloat fogColor[4] = { 0.3f, 0.3f, 0.3f, 1.0f };      // Fog Color
		glFogfv(GL_FOG_COLOR, fogColor);
	}
	else{
		glDisable(GL_FOG);
	}




	//Vale tin kamera
	glLoadIdentity();
	gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2],
		0, 0, 0,
		0, 1, 0);


	

	//Draw Skybox
	if (drawSkyboxFlag){
		drawSkybox();
	}

	// Draw Terrain
	drawTerrain();

	//Draw Sea
	if (drawSeaFlag){
		drawSea();
	}

	drawSun();





	//Add Light
	glEnable(GL_NORMALIZE);
	glEnable(GL_LINEAR_ATTENUATION);
	glShadeModel(GL_SMOOTH);
	float zeros[3] = { 0, 0, 0 };

	glMaterialfv(GL_FRONT, GL_EMISSION, zeros);
	GLfloat lightPos[] = { lightPosition[0], lightPosition[1], lightPosition[2], 1 };
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);


	glutSwapBuffers();
}





// Initialization Function
void myinit(void)
{

	//Load the skybox texture
	skyboxTexture = loadBMP_custom(eveningSkyboxFilename);
	elevationMap = loadBMP_custom(elevationFilename);
	seaTexture = loadBMP_custom(seaFilename);


	glClearColor(0, 0, 0, 0);


	int terrainDims = pow(2, terrainDepth) + 1;


	//Calculate the terrain
	point3  a = { -terrainSize, 0, terrainSize };
	point3 b = { terrainSize, 0, terrainSize };
	point3 c = { terrainSize, 0, -terrainSize };
	point3 d = { -terrainSize, 0, -terrainSize };
	calculateTerrain(a, b, c, d);




	// Specific spot effects
	// Cut off angle is 60 degrees
	glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 180.0f);
	// Fairly shiny spot
	glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 100.0f);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.000025);


	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	GLfloat args[] = { 1, 1, 1, 1 };

	glLightfv(GL_LIGHT0, GL_SPECULAR, args);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, args);
	args[0] = args[1] = args[2] = 0.1;
	glLightfv(GL_LIGHT0, GL_AMBIENT, args);

	//arxikh 8esh kameras
	cameraAngle = 45;
	cameraPosition[0] = cameraCenter[0] + cameraRadius * cos(cameraAngle * PI / 180);
	cameraPosition[1] = cameraCenter[1];
	cameraPosition[2] = cameraCenter[2] + cameraRadius * sin(cameraAngle * PI / 180);


	//arxikh 8esh fwtos
	lightAngle = 180;
	lightPosition[0] = lightCenter[0] + lightRadius * cos(lightAngle * PI / 180); //ousiastika x = -50
	lightPosition[1] = lightCenter[1]; // y = 0
	lightPosition[2] = lightCenter[2]; // z = 0





	//vale tin kamera
	glLoadIdentity();
	gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2],
		0, 0, 0,
		0, 1, 0);


	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fieldOfView, 1, 1, 1500);
	glMatrixMode(GL_MODELVIEW);


}


/* reshaped window */
void reshape(int width, int height) {
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fieldOfView, (GLfloat)width / (GLfloat)height, 1, 1500.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}



void main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(450, 0);
	glutCreateWindow("Fractal Terrain");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);

	int skyboxSubmenu = glutCreateMenu(RightClickMenu);
	glutAddMenuEntry("Day", 1);
	glutAddMenuEntry("Night", 2);
	glutAddMenuEntry("Evening", 3);
	glutAddMenuEntry("Interstellar", 4);
	glutAddMenuEntry("No skybox", 5);

	int qualitySubmenu = glutCreateMenu(RightClickMenu);
	glutAddMenuEntry("Low quality", 6);
	glutAddMenuEntry("Medium quality", 7);
	glutAddMenuEntry("High quality", 8);
	glutAddMenuEntry("Ultra High quality", 9);

	int roughnessSubmenu = glutCreateMenu(RightClickMenu);
	glutAddMenuEntry("Low Roughness", 10);
	glutAddMenuEntry("Medium Roughness", 11);
	glutAddMenuEntry("High Roughness", 12);


	glutCreateMenu(RightClickMenu);
	glutAddSubMenu("Skybox", skyboxSubmenu);
	glutAddSubMenu("Quality", qualitySubmenu);
	glutAddSubMenu("Roughness", roughnessSubmenu);
	glutAddMenuEntry("Show/Hide Sea", 13);
	glutAddMenuEntry("Show/Hide Wireframe", 14);
	glutAddMenuEntry("Show/Hide Fog", 15);
	glutAddMenuEntry("Show/Hide Heightmap", 16);
	glutAddMenuEntry("Exit", 17);
	glutAttachMenu(GLUT_RIGHT_BUTTON);



	glutKeyboardFunc(keyPressed);
	glutKeyboardUpFunc(keyReleased);
	glutSpecialFunc(specialKeyPressed);
	glutSpecialUpFunc(specialKeyReleased);
	glutIdleFunc(idle);

	myinit();



	glutMainLoop();
}