#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <queue>
using namespace std;

#define ImageW 1000
#define ImageH 700

float framebuffer[ImageH][ImageW][3];
GLFWwindow *pWindow;

struct color {
	float r, g, b;		// Color (R,G,B values)
	color(float rr, float gg, float bb) : r(rr), g(gg), b(bb) {}
	color() {}
	void set(float rr, float gg, float bb) {
		r = rr;
		g = gg;
		b = bb;
	}
	friend bool operator==(const color & lhs, const color & rhs) {
		return lhs.r == rhs.r && lhs.g == rhs.g && lhs.b == rhs.b;
	}
	friend bool operator!=(const color & lhs, const color & rhs) {
		return lhs.r != rhs.r && lhs.g != rhs.g && lhs.b != rhs.b;
	}
};

struct polyPoint {
	double x;
	double y;
	color col;
	polyPoint(double xx, double yy, color c) {
		if (xx <= 0) {
			x = 0;
		}
		if (yy <= 0) {
			y = 0;
		}
		x = xx;
		y = yy;
		col = c;
	}
	polyPoint(double xx, double yy) : x(xx), y(yy) {
		if (xx <= 0) {
			x = 0;
		}
		if (yy <= 0) {
			y = 0;
		}
		x = xx;
		y = yy;
	}
	polyPoint() {}
};

struct edge {
	double currentX;
	double slope;
	double maxY;
	polyPoint startpt;
	polyPoint endpt;
	edge(double cX, double s, double mY) : currentX(cX), slope(s), maxY(mY) {}
	edge(double cX, double s, double mY, polyPoint spt, polyPoint ept) : currentX(cX), slope(s), maxY(mY), startpt(spt), endpt(ept) {}
	void debug() {
		cout << "edge's currentX: " << currentX << " | slope: " << slope << " | maxY: " << maxY << endl;
	}
};

// GLOBAL VARIABLES
color curColor(1.0, 1.0, 1.0);
char drawmode = ' ';
vector<double> drawn_line_x;
vector<double> drawn_line_y;
double start_x, start_y, end_x, end_y;
vector<polyPoint> curPolygonPoints;
vector<edge> possibleEdges;

// Draws the scene
void drawit(void)
{
	glDrawPixels(ImageW, ImageH, GL_RGB, GL_FLOAT, framebuffer);
}

// Clears framebuffer to black
void clearFramebuffer()
{
	int i, j;

	for (i = 0; i<ImageH; i++) {
		for (j = 0; j<ImageW; j++) {
			framebuffer[i][j][0] = 0.0;
			framebuffer[i][j][1] = 0.0;
			framebuffer[i][j][2] = 0.0;
		}
	}
}

// Origin is upper left corner. 
void setFramebuffer(int x, int y, float R, float G, float B)
{
	// changes the origin from the lower-left corner to the upper-left corner
	y = ImageH - 1 - y;
	if (R <= 1.0)
		if (R >= 0.0)
			framebuffer[y][x][0] = R;
		else
			framebuffer[y][x][0] = 0.0;
	else
		framebuffer[y][x][0] = 1.0;
	if (G <= 1.0)
		if (G >= 0.0)
			framebuffer[y][x][1] = G;
		else
			framebuffer[y][x][1] = 0.0;
	else
		framebuffer[y][x][1] = 1.0;
	if (B <= 1.0)
		if (B >= 0.0)
			framebuffer[y][x][2] = B;
		else
			framebuffer[y][x][2] = 0.0;
	else
		framebuffer[y][x][2] = 1.0;
}

void display(void)
{
	// should not be necessary but some GPUs aren't displaying anything until a clear call.
	glClear(GL_COLOR_BUFFER_BIT);

	drawit();
	glFlush();
}

void dda(GLFWwindow* pwnd, double& start_x, double& start_y, double& end_x, double& end_y, color cur) {

	float step;
	double dx = (end_x - start_x);
	double dy = (end_y - start_y);
	if (abs(dx) >= abs(dy)) {
		step = abs(dx);
	}
	else {
		step = abs(dy);
	}
	dx = dx / step;
	dy = dy / step;
	double x = start_x;
	double y = start_y;
	int i = 1;
	while (i <= step) {
		setFramebuffer(x, y, cur.r, cur.g, cur.b);
		x += dx;
		y += dy;
		++i;
	}
}

vector<edge> merge(vector<edge> first, vector<edge> second) {
	vector<edge> merged;

	int i = 0, j = 0, k = 0;
	while (i < first.size() && j < second.size()) {
		if (first[i].currentX <= second[i].currentX) {
			merged.push_back(first[i]);
			++i;
		}
		else {
			merged.push_back(second[i]);
			++j;
		}
		++k;
	}

	while (i < first.size()) {
		merged.push_back(first[i]);
		++i;
		++k;
	}
	while (j < second.size()) {
		merged.push_back(second[j]);
		++j;
		++k;
	}

	return merged;
}

void mergesort(vector<edge>& edgeList) {

	if (edgeList.size() <= 1) {
		return;
	}

	vector<edge> firsthalf;
	vector<edge> secondhalf;
	for (int i = 0; i < edgeList.size() / 2; ++i) {
		firsthalf.push_back(edgeList[i]);
	}
	for (int i = edgeList.size() / 2; i < edgeList.size(); ++i) {
		secondhalf.push_back(edgeList[i]);
	}

	mergesort(firsthalf);
	mergesort(secondhalf);

	edgeList = merge(firsthalf, secondhalf);
}

void polyScanConversion(GLFWwindow* pwnd, vector<vector<edge>> activeEdgeTable, vector<edge> listOfEdges, color cur) {
	for (int i = 0; i < ImageH; ++i) {
		vector<edge> edgeToPush;
		for (int j = 0; j < ImageW; ++j) {
			for (int k = 0; k < listOfEdges.size(); ++k) {
				if (listOfEdges[k].currentX == j && ( listOfEdges[k].startpt.y == i || listOfEdges[k].endpt.y == i )) {
					vector<edge>::iterator it = edgeToPush.begin();
					if (listOfEdges[k].startpt.y == i) {
						listOfEdges[k].currentX = listOfEdges[k].startpt.x;
					}
					else {
						listOfEdges[k].currentX = listOfEdges[k].endpt.x;
					}
					edgeToPush.insert(it, listOfEdges[k]);
					listOfEdges.erase(listOfEdges.begin() + k);
				}
			}
		}
		
		activeEdgeTable[i] = edgeToPush;
	}

	vector<edge> activeEdgeList;

	for (int i = 0; i < activeEdgeTable.size(); ++i) {
		// add from AET to AEL
		for (int j = 0; j < activeEdgeTable[i].size(); ++j) {
			vector<edge>::iterator it = activeEdgeList.begin();
			activeEdgeList.insert(it, activeEdgeTable[i][j]);
		}
		cout << endl;
		// delete from AEL if the max value was passed
		for (int k = 0; k < activeEdgeList.size(); ++k) {
			if (activeEdgeList[k].maxY <= i) {
				activeEdgeList.erase(activeEdgeList.begin() + k);
			}
		}

		// sort by currentX
		mergesort(activeEdgeList);
		
		// fill pixels of pairs
		if (activeEdgeList.size() > 1) {
			for (int l = 0; l < activeEdgeList.size() - 1; ++l) {
				int amt_of_pixels = round(activeEdgeList[l + 1].currentX) - round(activeEdgeList[l].currentX);
				int px = (int)(round(activeEdgeList[l].currentX));
				while (amt_of_pixels > 0) {
					setFramebuffer(px, i, cur.r, cur.g, cur.b);
					++px;
					--amt_of_pixels;
				}
			}
			// increment current x
			for (int i = 0; i < activeEdgeList.size(); ++i) 
				activeEdgeList[i].currentX += activeEdgeList[i].slope;
		}
		// increment y, done by for loop
	}
}

bool sameColor(polyPoint pt1, polyPoint pt2) {
	if (ImageH - 1 - pt1.y < 0 || ImageH - 1 - pt2.y < 0) {
		cout << "out of bounds" << endl;
		return false;
	}

	return (framebuffer[ImageH-1-(int)pt1.y][(int)pt1.x][0] == framebuffer[ImageH - 1 - (int)pt2.y][(int)pt2.x][0] &&
		framebuffer[ImageH - 1 - (int)pt1.y][(int)pt1.x][1] == framebuffer[ImageH - 1 - (int)pt2.y][(int)pt2.x][1] &&
		framebuffer[ImageH - 1 - (int)pt1.y][(int)pt1.x][2] == framebuffer[ImageH - 1 - (int)pt2.y][(int)pt2.x][2]);
}

bool sameColor(polyPoint pt1, color cur) {
	if (ImageH - 1 - pt1.y < 0) {
		cout << "out of bounds" << endl;
		return false;
	}
	return (framebuffer[ImageH - 1 -(int)pt1.y][(int)pt1.x][0] == cur.r &&
		framebuffer[ImageH - 1 -(int)pt1.y][(int)pt1.x][1] == cur.g &&
		framebuffer[ImageH - 1- (int)pt1.y][(int)pt1.x][2] == cur.b);
}

void floodFill(GLFWwindow* pwnd, double ptX, double ptY, color cur, color initcolor) {
	if (cur == initcolor) {
		return;
	}
	if (!sameColor(polyPoint(ptX, ptY), initcolor)) {
		return;
	}
	queue<polyPoint> ptstack;
	polyPoint firstPt = polyPoint(ptX, ptY);
	setFramebuffer(firstPt.x, firstPt.y, cur.r, cur.g, cur.b);
	ptstack.push(firstPt);
	while (!ptstack.empty()) {
		polyPoint pt = ptstack.front();		
		ptstack.pop();
		polyPoint rightPt = polyPoint(pt.x + 1, pt.y);
		polyPoint leftPt = polyPoint(pt.x - 1, pt.y);
		polyPoint upPt = polyPoint(pt.x, pt.y + 1);
		polyPoint downPt = polyPoint(pt.x, pt.y - 1);
		if (sameColor(rightPt, initcolor)) {
			setFramebuffer(rightPt.x, rightPt.y, cur.r, cur.g, cur.b);
			ptstack.push(rightPt);
		}
		if (sameColor(leftPt, initcolor)) {
			setFramebuffer(leftPt.x, leftPt.y, cur.r, cur.g, cur.b);
			ptstack.push(leftPt);
		}
		if (sameColor(downPt, initcolor)) {
			setFramebuffer(downPt.x, downPt.y, cur.r, cur.g, cur.b);
			ptstack.push(downPt);
		}
		if (sameColor(upPt, initcolor)) {
			setFramebuffer(upPt.x, upPt.y, cur.r, cur.g, cur.b);
			ptstack.push(upPt);
		}
	}
	return;
}

void boundaryFill(GLFWwindow* pwnd, double ptX, double ptY, color cur, color boundcolor) {
	if (!sameColor(polyPoint(ptX, ptY), boundcolor) && !sameColor(polyPoint(ptX, ptY), cur)) {

		setFramebuffer(ptX, ptY, cur.r, cur.g, cur.b);
		boundaryFill(pwnd, ptX + 1, ptY, cur, boundcolor);
		boundaryFill(pwnd, ptX - 1, ptY, cur, boundcolor);
		boundaryFill(pwnd, ptX, ptY + 1, cur, boundcolor);
		boundaryFill(pwnd, ptX, ptY - 1, cur, boundcolor);
	}
	return;
}

void freeHand(GLFWwindow* pwnd, double xpos, double ypos) {
	if (drawn_line_x.size() > 1) {
		for (int i = 0; i < drawn_line_x.size() - 1; ++i) {
			dda(pwnd, drawn_line_x[i], drawn_line_y[i], drawn_line_x[i + 1], drawn_line_y[i + 1], curColor);
		}
	}
	drawn_line_x.push_back(xpos);
	drawn_line_y.push_back(ypos);
}

double edgeSlope(double x1, double y1, double x2, double y2) {
	double initial = (y2 - y1) / (x2 - x1);
	return 1 / initial;
}

double findMaxY(double y1, double y2) {
	if (y1 > y2) return y1;
	else return y2;
}

void drawLineCallback(GLFWwindow* pwnd, double xpos, double ypos) {
	drawn_line_x.push_back(xpos);
	drawn_line_y.push_back(ypos);
	for (int i = 0; i < drawn_line_x.size(); ++i) {
		dda(pwnd, start_x, start_y, drawn_line_x[i], drawn_line_y[i], color(0.0,0.0,0.0));			
	}
	dda(pwnd, start_x, start_y, xpos, ypos, curColor);
}

// Mouse action callback function.
// button can be GLFW_MOUSE_BUTTON_LEFT or GLFW_MOUSE_BUTTON_RIGHT
// action can be GLFW_PRESS or GLFW_RELEASE
void MouseButtonCallbackFunc(GLFWwindow* pwnd, int button, int action, int mods)
{

	switch (drawmode)
	{
	case 'L':
	{
		// line drawing mode
		// TODO: better solution than make x and y variables global
		if (action == GLFW_PRESS) {
			glfwGetCursorPos(pwnd, &start_x, &start_y);
			glfwSetCursorPosCallback(pwnd, drawLineCallback);
		}
		else if (action == GLFW_RELEASE) {
			//glfwGetCursorPos(pwnd, &end_x, &end_y);
			//dda(pwnd, start_x, start_y, end_x, end_y, curColor);
			glfwSetCursorPosCallback(pwnd, NULL);
			drawn_line_x = {};
			drawn_line_y = {};
		}
		break;
	}
	case 'P':
	{
		// polygon drawing mode
		vector<vector<edge>> activeEdgeTable(ImageH);
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			double newx, newy;
			glfwGetCursorPos(pwnd, &newx, &newy);
			curPolygonPoints.push_back(polyPoint(newx, newy));
			cout << "Just added x of " << newx << " and y of " << newy << endl;
		}
		else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
			double finalx, finaly;
			glfwGetCursorPos(pwnd, &finalx, &finaly);
			curPolygonPoints.push_back(polyPoint(finalx, finaly));
			for (int i = 0; i < curPolygonPoints.size() - 1; ++i) {
				dda(pwnd, curPolygonPoints[i].x, curPolygonPoints[i].y, curPolygonPoints[i + 1].x, curPolygonPoints[i + 1].y, curColor);
				edge newEdge = edge(
										curPolygonPoints[i].x,
										edgeSlope(curPolygonPoints[i].x, curPolygonPoints[i].y,
											curPolygonPoints[i + 1].x, curPolygonPoints[i + 1].y),
										findMaxY(curPolygonPoints[i].y, curPolygonPoints[i + 1].y),
										polyPoint(curPolygonPoints[i].x, curPolygonPoints[i].y),
										polyPoint(curPolygonPoints[i + 1].x, curPolygonPoints[i + 1].y)
									);
				newEdge.debug();
				possibleEdges.push_back(newEdge);
			}
			int last = curPolygonPoints.size() - 1;
			dda(pwnd, curPolygonPoints[0].x, curPolygonPoints[0].y, curPolygonPoints[last].x, curPolygonPoints[last].y, curColor);
			edge lastEdge = edge(
										curPolygonPoints[0].x,
										edgeSlope(curPolygonPoints[0].x, curPolygonPoints[0].y,
											curPolygonPoints[last].x, curPolygonPoints[last].y),
										findMaxY(curPolygonPoints[0].y, curPolygonPoints[last].y),
										polyPoint(curPolygonPoints[0].x, curPolygonPoints[0].y),
										polyPoint(curPolygonPoints[last].x, curPolygonPoints[last].y)
								);
			lastEdge.debug();
			possibleEdges.push_back(lastEdge);
			cout << possibleEdges.size() << endl;
			polyScanConversion(pwnd, activeEdgeTable, possibleEdges, curColor);
			possibleEdges = {};
			activeEdgeTable.clear();
			curPolygonPoints = {};

		}
		break;
	}
	case 'F':
	{
		// free drawing mode
		if (action == GLFW_PRESS) {
			double newx, newy;
			glfwGetCursorPos(pwnd, &newx, &newy);
			drawn_line_x.push_back(newx);
			drawn_line_y.push_back(newy);
			glfwSetCursorPosCallback(pwnd, freeHand);
		}
		if (action == GLFW_RELEASE) {
			glfwSetCursorPosCallback(pwnd, NULL);
			drawn_line_x = {};
			drawn_line_y = {};
		}
		break;
	}
	case 'Q':
	{
		// boundary fill mode
		double newx, newy;
		glfwGetCursorPos(pwnd, &newx, &newy);
		boundaryFill(pwnd, newx, newy, curColor, curColor);
		break;
	}
	case 'O':
	{
		// flood fill mode
		double newx, newy;
		glfwGetCursorPos(pwnd, &newx, &newy);
		cout << newx << " " << newy << endl;
		color thisCol = color(framebuffer[ImageH-1-(int)newy][(int)newx][0], framebuffer[ImageH-1-(int)newy][(int)newx][1],
			framebuffer[ImageH-1-(int)newy][(int)newx][2]);
		cout << "color " << thisCol.r << ", " << thisCol.g << ", " << thisCol.b << endl;
		floodFill(pwnd, newx, newy, curColor, thisCol);
	}
	}
}

// Keyboard action callback function.
// key is the action button. Use capital letters to compare. For example,
// if(key == 'A') 
// {
//		//.............
// }
// or
// switch(key)
// {
//		case 'A':
//			//....
//			break;
//		case 'B':
//			//...
//			break;
//			//....
// }
//
// action can be GLFW_PRESS or GLFW_RELEASE
void KeyCallbackFunc(GLFWwindow* pwnd, int key, int scancode, int action, int mode)
{
	if (action == GLFW_PRESS)
		std::cout << "Key " << (char)key << " is pressed." << std::endl;
	else if(action == GLFW_RELEASE)
		std::cout << "Key " << (char)key << " is released." << std::endl;
	
	switch ((char)key)
	{
		case 'C':
		{
			clearFramebuffer();
			break;
		}
		case '0': 
		{
			curColor.set(0.0, 0.0, 0.0);
			break;
		}
		case '1':
		{
			curColor.set(1.0, 0.0, 0.0);
			break;
		}
		case '2':
		{
			curColor.set(0.0, 1.0, 0.0);
			break;
		}
		case '3':
		{
			curColor.set(1.0, 1.0, 0.0);
			break;
		}
		case '4':
		{
			curColor.set(0.0, 0.0, 1.0);
			break;
		}
		case '5':
		{
			curColor.set(1.0, 0.0, 1.0);
			break;
		}
		case '6':
		{
			curColor.set(0.0, 1.0, 1.0);
			break;
		}
		case '7':
		{
			curColor.set(1.0, 1.0, 1.0);
			break;
		}
		case 'R':
		{
			curColor.set(curColor.r + 0.05f, curColor.g, curColor.b);
			break;
		}
		case 'G':
		{
			curColor.set(curColor.r, curColor.g + 0.05f, curColor.b);
			break;
		}
		case 'B':
		{
			curColor.set(curColor.r, curColor.g, curColor.b + 0.05f);
			break;
		}
		case 'L':
		{
			drawmode = 'L';
			break;
		}
		case 'P':
		{
			drawmode = 'P';
			break;
		}
		case 'F':
		{
			drawmode = 'F';
			break;
		}
		case 'Q':
		{
			drawmode = 'Q';
			break;
		}
		case 'O':
		{
			drawmode = 'O';
			break;
		}
	}
}

void init()
{
	glfwInit();
	glfwSetTime(0.0);
	if (pWindow)
	{
		glfwDestroyWindow(pWindow);
	}
	pWindow = glfwCreateWindow(ImageW, ImageH, "Assignment 1 - Aksel Taylan", NULL, NULL);
	glfwMakeContextCurrent(pWindow);
	glfwSetMouseButtonCallback(pWindow, MouseButtonCallbackFunc);
	glfwSetKeyCallback(pWindow, KeyCallbackFunc);
	glewExperimental = true;
	glewInit();
	glViewport(0, 0, ImageW, ImageH);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	clearFramebuffer();
}


int main()
{
	init();
	while (!glfwWindowShouldClose(pWindow)) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		display();
		glfwSwapBuffers(pWindow);
		glfwPollEvents();
	}
	return 0;
}