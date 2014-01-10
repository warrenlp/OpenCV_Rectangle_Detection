// The "Square Detector" program.
// It loads several images sequentially and tries to find squares in
// each image

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <iostream>
#include <math.h>
#include <string.h>

using namespace cv;
using namespace std;

static void help()
{
    cout <<
    "\nA program using pyramid scaling, Canny, contours, contour simpification and\n"
    "memory storage (it's got it all folks) to find\n"
    "squares in a list of images pic1-6.png\n"
    "Returns sequence of squares detected on the image.\n"
    "the sequence is stored in the specified memory storage\n"
    "Call:\n"
    "./squares\n"
    "Using OpenCV version %s\n" << CV_VERSION << "\n" << endl;
}


int thresh = 50, N = 11;
const char* wndname = "Square Detection Demo";

struct Center
{
    Point2d location;
    double radius;
    // double confidence;
};

// helper function:
// finds a cosine of angle between vectors
// from pt0->pt1 and from pt0->pt2
static double angle( Point pt1, Point pt2, Point pt0 )
{
    double dx1 = pt1.x - pt0.x;
    double dy1 = pt1.y - pt0.y;
    double dx2 = pt2.x - pt0.x;
    double dy2 = pt2.y - pt0.y;
    return (dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10);
}

// filters out squares found based on color, position, and size
static void filterSquares(const vector<vector<Point> > colorSquares[], vector<vector<Point> > &squares)
{
	const int numOfColorPlanes = 3;
	const int maxsize = 100000;
	const int minDistBetweenSquares = 1000;
	vector<vector<Point> > filteredSquares[3];

	// Eliminate squares that are greater than maximum allowable size (min size has already been filtered out)
	for (int i=0; i < numOfColorPlanes; ++i)
	{
		vector<vector<Point> >::const_iterator citer = colorSquares[i].begin();
		for (; citer != colorSquares[i].end(); ++citer)
		{
			if (fabs(contourArea(Mat(*citer))) < maxsize)
			{
				filteredSquares[i].push_back(*citer);
			}
		}

	}

	// Sort squares based on their relative centers as well as size
	vector<vector<Center> > centers;
	vector<vector<Center> > newCenters;
	vector<Center> curCenters;

	for (size_t c=0; c<numOfColorPlanes; ++c)
	{
		// current filtered square color plane
		const vector<vector<Point> > &curFSCP = filteredSquares[c];

		vector<vector<Point> >::const_iterator curSquareIter = curFSCP.begin();
		for (; curSquareIter != curFSCP.end(); ++curSquareIter)
		{
			if (fabs(contourArea(Mat(*citer))) < maxsize)
			{
				filteredSquares[i].push_back(*citer);
			}
		}

		for (size_t i=0; i<curRects.size(); ++i)
		{
			bool isNew = true;
			Center curCenter;
			curCenter.location = Point2f((curFSCP),());
			for (size_t j = 0; j<centers.size() && !isNew; ++j)
			{
				double dist = norm(centers[j][ centers[j].size() / 2 ].location - curCenters[i].location);
				isNew = dist >= minDistBetweenSquares && dist >= centers[j][ centers[j].size() / 2 ].radius && dist >= curCenters[i].radius;
				if (!isNew)
				{
					centers[j].push_back(curCenters[i]);

					size_t k = centers[j].size() - 1;
					while( k > 0 && centers[j][k].radius < centers[j][k-1].radius )
					{
						centers[j][k] = centers[j][k-1];
						k--;
					}
					centers[j][k] = curCenters[i];
				}
			}
			if (isNew)
			{
				newCenters.push_back(vector<Center> (1, curCenters[i]));
				//centers.push_back(vector<Center> (1, curCenters[i]));
			}
		}

	}

	// Eliminate squares with the same relative center within color plane
//	for (int i=3; i < numOfColorPlanes; ++i)
//	{
//
//		for ()
//		{
//
//
//			// Adjust to have average center and size
//
//		}
//	}


	// Eliminate based on color


	// Eliminate based on relative center between colors

	for (int i=0; i < numOfColorPlanes; ++i)
	{
		std::copy(filteredSquares[i].begin(), filteredSquares[i].end(), std::back_inserter(squares));
	}
}

// returns sequence of squares detected on the image.
// the sequence is stored in the specified memory storage
static void findSquares( const Mat& image, vector<vector<Point> >& squares )
{
    squares.clear();

    Mat pyr, timg, gray0(image.size(), CV_8U), gray;

    // down-scale and upscale the image to filter out the noise
    pyrDown(image, pyr, Size(image.cols/2, image.rows/2));
    pyrUp(pyr, timg, image.size());
    vector<vector<Point> > contours;
    vector<vector<Point> > colorSquares[3];

    // find squares in every color plane of the image
    for( int c = 0; c < 3; c++ )
    {
    	colorSquares[c].clear();

        int ch[] = {c, 0};
        mixChannels(&timg, 1, &gray0, 1, ch, 1);

        // try several threshold levels
        for( int l = 0; l < N; l++ )
        {
            // hack: use Canny instead of zero threshold level.
            // Canny helps to catch squares with gradient shading
            if( l == 0 )
            {
                // apply Canny. Take the upper threshold from slider
                // and set the lower to 0 (which forces edges merging)
                Canny(gray0, gray, 0, thresh, 5);
                // dilate canny output to remove potential
                // holes between edge segments
                dilate(gray, gray, Mat(), Point(-1,-1));
            }
            else
            {
                // apply threshold if l!=0:
                //     tgray(x,y) = gray(x,y) < (l+1)*255/N ? 255 : 0
                gray = gray0 >= (l+1)*255/N;
            }

            // find contours and store them all as a list
            findContours(gray, contours, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);

            vector<Point> approx;

            // test each contour
            for( size_t i = 0; i < contours.size(); i++ )
            {
                // approximate contour with accuracy proportional
                // to the contour perimeter
                approxPolyDP(Mat(contours[i]), approx, arcLength(Mat(contours[i]), true)*0.02, true);

                // square contours should have 4 vertices after approximation
                // relatively large area (to filter out noisy contours)
                // and be convex.
                // Note: absolute value of an area is used because
                // area may be positive or negative - in accordance with the
                // contour orientation
                if( approx.size() == 4 &&
                    fabs(contourArea(Mat(approx))) > 1000 &&
                    isContourConvex(Mat(approx)) )
                {
                    double maxCosine = 0;

                    for( int j = 2; j < 5; j++ )
                    {
                        // find the maximum cosine of the angle between joint edges
                        double cosine = fabs(angle(approx[j%4], approx[j-2], approx[j-1]));
                        maxCosine = MAX(maxCosine, cosine);
                    }

                    // if cosines of all angles are small
                    // (all angles are ~90 degree) then write quandrange
                    // vertices to resultant sequence
                    if( maxCosine < 0.3 )
                        colorSquares[c].push_back(approx);
                }
            }
        }
    }

    filterSquares(colorSquares, squares);
}


// the function draws all the squares in the image
static void drawSquares( Mat& image, const vector<vector<Point> >& squares )
{
    for( size_t i = 0; i < squares.size(); i++ )
    {
        const Point* p = &squares[i][0];
        int n = (int)squares[i].size();
        polylines(image, &p, &n, 1, true, Scalar(0,255,0), 3, CV_AA);
    }

    imshow(wndname, image);
}


int main(int /*argc*/, char** /*argv*/)
{
//    static const char* names[] = { "pic1.png", "pic2.png", "pic3.png",
//        "pic4.png", "pic5.png", "pic6.png", 0 };
	static const char* names[] = { "PaperRectangle.jpg", "PaperRectangle2.jpg", "PaperRectangle_light.jpg", 0 };
    help();
    namedWindow( wndname, 1 );
    vector<vector<Point> > squares;

    for( int i = 0; names[i] != 0; i++ )
    {
        Mat image = imread(names[i], 1);
        if( image.empty() )
        {
            cout << "Couldn't load " << names[i] << endl;
            continue;
        }

        findSquares(image, squares);
        drawSquares(image, squares);

        int c = waitKey();
        if( (char)c == 27 )
            break;
    }

    return 0;
}
