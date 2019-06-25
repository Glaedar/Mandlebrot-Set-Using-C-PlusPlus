// Architectures and Performance: Mandelbrot
// Andrew Peters Coursework 1502220
//Data Structures and Algorithms 2


#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <thread>
#include <vector>
#include <functional>
#include <mutex>

// Import things we need from the standard library
using namespace std;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::complex;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::thread;
using std::mem_fun;
using std::unique_lock;
using std::condition_variable;

mutex task_mutex;
bool task_ready = false;
condition_variable task_cv;


// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_clock;


// The size of the image to generate.
const int WIDTH = 2048;
const int HEIGHT = 1152;

//the begining and the end of the first slice
int yStart = 0;
int yEnd = HEIGHT / 16;

// The number of times to iterate before we assume that a point isn't in the
// Mandelbrot set.
// (You may need to turn this up if you zoom further into the set.)
const int MAX_ITERATIONS = 500;

// The image data.
// Each pixel is represented as 0xRRGGBB.
uint32_t image[HEIGHT][WIDTH];


// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char *filename)
{
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *) header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char *) pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}

//creates a virtual class
class Task
{
public:
	virtual ~Task()
	{
	}
	virtual void run() = 0;
};

class MessageTask : public Task
{
public:

	// initialise the variables
	double xleft;
	double xright;
	double ytop;
	double ybottom;
	int start;
	int end;

	MessageTask(double left, double right, double top, double bottom, int yStart, int yEnd)
	{
		xleft = left;
		xright = right;
		ytop = top;
		ybottom = bottom;
		start = yStart;
		end = yEnd;
	}
	void run();

private:
};


class Farm
{
public:
	void add_task(Task *task);
	void run();
	void push_task();

private:
	void threadWorker();
	std::queue<Task *> task_queue;

	int numThreads = 4;
	int thread_divide = (numThreads * numThreads);

};


void MessageTask::run()
{
	//run the mandlebrot slice
	for (int y = start; y < end; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			complex<double> c(xleft + (x * (xright - xleft) / WIDTH),
				ytop + (y * (ybottom - ytop) / HEIGHT));

			// Start off z at (0, 0).
			complex<double> z(0.0, 0.0);

			// Iterate z = z^2 + c until z moves more than 2 units
			// away from (0, 0), or we've iterated too many times.
			int iterations = 0;
			while (abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = (z * z) + c;

				++iterations;
			}

			if (iterations == MAX_ITERATIONS)
			{
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = 0x000000; // black
			}
			else
			{
				// z escaped within less than MAX_ITERATIONS
				// iterations. This point isn't in the set.
				image[y][x] = (iterations * 5 << 6) | (iterations * 5 << 12) | (iterations * 5 << 7); // white
			}
		}
	}

}
void Farm::add_task(Task *task)
{

	//add a task and notify a sleeping thread
	unique_lock<mutex> lock(task_mutex);
	task_queue.push(task);
	task_ready = true;
	task_cv.notify_one();


	
};

void Farm::run()
{

	std::vector<thread *> farmthreads;

	//create multiple threads using the 2 thread functions
	for (int i = 0; i < numThreads; i++)
	{
		farmthreads.push_back(new thread(mem_fun(&Farm::push_task), this));

		farmthreads.push_back(new thread(mem_fun(&Farm::threadWorker), this));

	}
	for (thread *t: farmthreads)
	{
		t->join();
		delete t;
		
	}
}
void Farm ::threadWorker()
{
	//checks to see if the task queue is empty or not and then runs the front of the queue
	
	while (!task_queue.empty())
	{
		Task* nextTask;
		{
			unique_lock<mutex> lock(task_mutex);
			
			while (!task_ready)
			{
				task_cv.wait(lock);
			}
			nextTask = task_queue.front();
			task_queue.pop();

		}
		nextTask->run();
		delete nextTask;
	}
	
}

void Farm::push_task()
{

	//puts the task into the task queue
	for (int i = 0; i < numThreads; i++)
	{
		add_task(new MessageTask(-2.0, 1.0, 1.125, -1.125, yStart, yEnd));

		yStart = yEnd;
		yEnd = (HEIGHT / thread_divide) + yEnd;

	}
}
int main(int argc, char *argv[])
{
	//create an excel file to hold the data of times when it completes
	std::ofstream computeTime;
	computeTime.open("2048x1124.csv", std::ios_base::app);

	// Start timing
	the_clock::time_point start = the_clock::now();
	

	std::vector<thread *> taskthreads;
	std::vector<thread*> threadNum;

	Farm farm;

	farm.run();

	// Stop timing
	the_clock::time_point end = the_clock::now();

	// Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Computing the Mandelbrot set took " << time_taken << " ms." << endl;
	
	computeTime << WIDTH << "," << HEIGHT<< "," << time_taken << endl;
	write_tga("output.tga");

	computeTime.close();	

	return 0;
}
