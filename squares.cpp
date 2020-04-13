// The "Square Detector" program.
// It loads several images sequentially and tries to find squares in
// each image

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <iostream>
#include <math.h>
#include <string.h>
#include <numeric>

using namespace cv;
using namespace std;

static void help()
{
    cout << "\nA program using pyramid scaling, Canny, contours, contour simpification and\n"
            "memory storage (it's got it all folks) to find\n"
            "squares in a list of images pic1-6.png\n"
            "Returns sequence of squares detected on the image.\n"
            "the sequence is stored in the specified memory storage\n"
            "Call:\n"
            "./squares\n"
            "Using OpenCV version %s\n"
         << CV_VERSION << "\n"
         << endl;
}

int thresh = 50, N = 11;
const char *wndname = "Square Detection Demo";

struct Center
{
    Point2f location;
    float radius;
    // double confidence;
};

// helper function:
// finds a cosine of angle between vectors
// from pt0->pt1 and from pt0->pt2
static double angle(Point pt1, Point pt2, Point pt0)
{
    double dx1 = pt1.x - pt0.x;
    double dy1 = pt1.y - pt0.y;
    double dx2 = pt2.x - pt0.x;
    double dy2 = pt2.y - pt0.y;
    return (dx1 * dx2 + dy1 * dy2) / sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2) + 1e-10);
}

//*****************************************************************************
/*!
 *  \brief  FilterByMaxSize
 *
 *  \param  colorSquares, filteredBySizeSquares, maxsize
 *
 *  \version
 *      - W Parsons   01/12/2014
 *        Initial Version
 *
 *****************************************************************************/
static void FilterByMaxSize(const vector<vector<Point>> &colorSquares, vector<vector<Point>> &filteredBySizeSquares, const int maxsize = 100000)
{
    // Eliminate squares that are greater than maximum allowable size (min size has already been filtered out)
    vector<vector<Point>>::const_iterator citer = colorSquares.begin();
    for (; citer != colorSquares.end(); ++citer)
    {
        if (fabs(contourArea(Mat(*citer))) < maxsize)
        {
            filteredBySizeSquares.push_back(*citer);
        }
    }
}

//*****************************************************************************
/*!
 *  \brief  FilterByMaxSize
 *
 *  \version
 *      - W Parsons   01/12/2014
 *        Initial Version
 *
 *****************************************************************************/
static void FilterByMaxSize(const vector<vector<Point>> colorSquares[], vector<vector<Point>> filteredBySizeSquares[], const int maxsize = 100000)
{
    const int numOfColorPlanes = 3;
    // Eliminate squares that are greater than maximum allowable size (min size has already been filtered out)
    for (int c = 0; c < numOfColorPlanes; ++c)
    {
        FilterByMaxSize(colorSquares[c], filteredBySizeSquares[c], maxsize);
    }
}

//*****************************************************************************
/*!
 *  \brief  SortByCenters
 *
 *  \param  colorSquares, sortedSquares, minDist
 *
 *  \version
 *      - W Parsons   01/12/2014
 *        Initial Version
 *
 *****************************************************************************/
static void SortByCenters(const vector<vector<Point>> &colorSquares, vector<vector<vector<Point>>> &sortedSquares,
                          const int minDist = 50)
{
    for (size_t i = 0; i < colorSquares.size(); ++i)
    {
        bool isNew = true;
        Center curCenter;
        minEnclosingCircle(Mat(colorSquares[i]), curCenter.location, curCenter.radius);

        for (size_t j = 0; j < sortedSquares.size() && isNew; ++j)
        {
            Center sortedCenter;

            // For this algorithm we pick the sorted square in the middle of the array to check against since they are sorted by radius size
            minEnclosingCircle(Mat(sortedSquares[j][sortedSquares[j].size() / 2]), sortedCenter.location, sortedCenter.radius);

            double dist = norm(sortedCenter.location - curCenter.location);
            isNew = dist >= minDist && dist >= sortedCenter.radius && dist >= curCenter.radius;
            if (!isNew)
            {
                // Determine where this radius fits in the group
                size_t k = sortedSquares[j].size() - 1;
                minEnclosingCircle(Mat(sortedSquares[j][k]), sortedCenter.location, sortedCenter.radius);
                while (k > 0 && curCenter.radius < sortedCenter.radius)
                {
                    k--;
                    minEnclosingCircle(Mat(sortedSquares[j][k]), sortedCenter.location, sortedCenter.radius);
                }
                if (curCenter.radius > sortedCenter.radius)
                {
                    ++k;
                }
                sortedSquares[j].insert(sortedSquares[j].begin() + k, colorSquares[i]);
            }
        }
        if (isNew)
        {
            // Start a new group of squares
            sortedSquares.push_back(vector<vector<Point>>(1, colorSquares[i]));
        }
    }
}

//*****************************************************************************
/*!
 *  \brief  SortByCenters
 *
 *  \version
 *      - W Parsons   01/12/2014
 *        Initial Version
 *
 *****************************************************************************/
static void SortByCenters(const vector<vector<Point>> colorSquares[], vector<vector<vector<Point>>> sortededByCentersSquares[],
                          const int minDist = 50)
{
    const int numOfColorPlanes = 3;
    // Eliminate squares that are greater than maximum allowable size (min size has already been filtered out)
    for (int c = 0; c < numOfColorPlanes; ++c)
    {
        SortByCenters(colorSquares[c], sortededByCentersSquares[c], minDist);
    }
}

//*****************************************************************************
/*!
 *  \brief  DumpSortedCenterData
 *
 *  \version
 *      - W Parsons   01/12/2014
 *        Initial Version
 *
 *****************************************************************************/
static void DumpSortedCenterData(const vector<vector<vector<Point>>> &sortedByCenterSquares)
{
    vector<vector<vector<Point>>>::const_iterator cCenterDataIter = sortedByCenterSquares.begin();
    for (; cCenterDataIter != sortedByCenterSquares.end(); ++cCenterDataIter)
    {
        vector<vector<Point>>::const_iterator cCenterIter = cCenterDataIter->begin();
        for (; cCenterIter != cCenterDataIter->end(); ++cCenterIter)
        {
            Center center;
            minEnclosingCircle(Mat(*cCenterIter), center.location, center.radius);
            cout << "L: " << center.location << " R: " << center.radius << "\t";
        }
        cout << endl;
    }

    cout << endl;
}

//*****************************************************************************
/*!
 *  \brief  DumpSortedCenterData
 *
 *  \version
 *      - W Parsons   01/12/2014
 *        Initial Version
 *
 *****************************************************************************/
static void DumpSortedCenterData(const vector<vector<vector<Point>>> sortedByCenterSquares[])
{
    const int numOfColorPlanes = 3;
    // Eliminate squares that are greater than maximum allowable size (min size has already been filtered out)
    for (int c = 0; c < numOfColorPlanes; ++c)
    {
        cout << "Dumping Data from layer \"" << c << "\" sortedByCenterSquares:" << endl;
        DumpSortedCenterData(sortedByCenterSquares[c]);
    }
}

//*****************************************************************************
/*!
 *  \brief  ConsolidateSquares
 *
 *  \version
 *      - W Parsons   01/13/2014
 *        Initial Version
 *
 *****************************************************************************/
static void ConsolidateSquares(const vector<vector<vector<Point>>> sortedSquares[], vector<vector<Point>> &consolidatedSquares)
{
    // Make sure we don't have any residual data in this vector or we could end up with repeats
    consolidatedSquares.clear();

    const int numOfColorPlanes = 3;
    vector<vector<vector<Point>>> sortedByCenterSquares;
    const int minDist = 50;

    // For this algorithm we simply select the square located in the middle of the vector since they are sorted by radius size

    // Pull center from each plane into one vector
    for (int c = 0; c < numOfColorPlanes; ++c)
    {
        vector<vector<vector<Point>>>::const_iterator cSortedIter = sortedSquares[c].begin();
        for (; cSortedIter != sortedSquares[c].end(); ++cSortedIter)
        {
            size_t middle = cSortedIter->size() / 2;
            vector<Point> middleSquare = cSortedIter->at(middle);
            consolidatedSquares.push_back(middleSquare);
        }
    }

    // Reduce to one unique square
    //	Color information is gone by this point
    SortByCenters(consolidatedSquares, sortedByCenterSquares, minDist);

    DumpSortedCenterData(sortedByCenterSquares);

    consolidatedSquares.clear();

    vector<vector<vector<Point>>>::const_iterator cSortedIter = sortedByCenterSquares.begin();
    for (; cSortedIter != sortedByCenterSquares.end(); ++cSortedIter)
    {
        size_t middle = cSortedIter->size() / 2;
        vector<Point> middleSquare = cSortedIter->at(middle);
        consolidatedSquares.push_back(middleSquare);
    }
}

//*****************************************************************************
/*!
 *  \brief  FilterByBGR
 *
 *  \version
 *      - W Parsons   01/13/2014
 *        Initial Version
 *
 *****************************************************************************/
static void FilterByBGR(const Mat &image, const vector<vector<Point>> &sortedSquares,
                        vector<vector<Point>> &sortedByRGBSquares,
                        const vector<Scalar> &colorRange)
{
    // Make sure that we don't have any residual data in the returned vector
    sortedByRGBSquares.clear();

    Mat upperMat;
    Mat lowerMat;
    Mat resultMat;

    bool isInRange(true);

    vector<vector<Point>>::const_iterator cSortedSquaresIter = sortedSquares.begin();
    for (; cSortedSquaresIter != sortedSquares.end(); ++cSortedSquaresIter)
    {
        Rect rect = boundingRect(*cSortedSquaresIter);

        // Reduce rect to 25% (50% per side) with center as anchor point
        Point2i topLeft(rect.tl());
        Point2i botRight(rect.br());
        Point2i delta = botRight - topLeft;
        topLeft += (delta * 0.25);
        botRight -= (delta * 0.25);

        rect = Rect(topLeft, botRight);

        // Define the Region of Interest (ROI)
        Mat imageRect = image(rect);

        lowerMat = Mat(imageRect.size(), CV_8UC3, colorRange[0]);
        upperMat = Mat(imageRect.size(), CV_8UC3, colorRange[1]);
        resultMat = Mat::zeros(imageRect.size(), CV_8U);

        inRange(imageRect, lowerMat, upperMat, resultMat);

        vector<Mat> planes;
        split(imageRect, planes);
        for (int c = 0; c < 3; ++c)
        {
            resultMat.copyTo(planes[c]);
        }
        merge(planes, imageRect);

        // For this algorithm, we make sure that we have 100% containment within the color range
        // 	For this we use a pixel-wise logical-AND operation
        int nr = resultMat.rows;
        int nc = resultMat.cols; // * resultMat.channels(); // Should be only one channel, but just in case...
        for (int i = 0; i < nr && isInRange; ++i)
        {
            uchar *data = resultMat.ptr<uchar>(i);
            for (int j = 0; j < nc && isInRange; ++j)
            {
                if (data[j] != 255)
                {
                    isInRange = false;
                }
            }
        }

        if (isInRange)
        {
            sortedByRGBSquares.push_back(*cSortedSquaresIter);
        }
        else
        {
            isInRange = true;
        }
    }
}

// filters out squares found based on color, position, and size
static void FilterSquares(const Mat &image, vector<vector<Point>> colorSquares[], vector<vector<Point>> &squares)
{
    squares.clear();

    const int maxsize = 100000;
    const int minDstBtnCtrs = 50;
    vector<vector<Point>> filteredBySizeSquares[3];
    vector<vector<vector<Point>>> sortedByCenterSquares[3];
    vector<vector<Point>> consolidatedSquares;
    vector<vector<Point>> sortedByRGBSquares;

    // Eliminate squares that are greater than maximum allowable size (min size has already been filtered out)
    FilterByMaxSize(colorSquares, filteredBySizeSquares, maxsize);

    // Sorts Squares into groups which share the same approximate center
    //	This is calculated based on the containing circle of the contours as well as the distance between centers
    SortByCenters(filteredBySizeSquares, sortedByCenterSquares, minDstBtnCtrs);

    DumpSortedCenterData(sortedByCenterSquares);

    // Reduces the squares to a single, unique square across all color planes
    ConsolidateSquares(sortedByCenterSquares, consolidatedSquares);

    // This is a very permissive range being the whole upper-half of the gray-scale spectrum
    //		Later on the algorithm will reject any image with any pixels out of this range, so we need it to be permissive
    vector<Scalar> colorRange;
    colorRange.push_back(Scalar(125, 125, 125));
    colorRange.push_back(Scalar(255, 255, 255));

    // Removes squares that do not fit within the prescribed color range
    FilterByBGR(image, consolidatedSquares, sortedByRGBSquares, colorRange);

    //TODO: Add FilterByHSL()???   // This might come in useful as resistance to different lighting conditions

    for (size_t c = 0; c < 3; ++c)
    {
        // Make sure we don't need the old, unfiltered squares data here any more
        colorSquares[c].clear();
        std::copy(filteredBySizeSquares[c].begin(), filteredBySizeSquares[c].end(), std::back_inserter(colorSquares[c]));
    }

    std::copy(sortedByRGBSquares.begin(), sortedByRGBSquares.end(), std::back_inserter(squares));
}

// returns sequence of squares detected on the image.
// the sequence is stored in the specified memory storage
static void FindSquares(const Mat &image, vector<vector<Point>> &squares, vector<vector<Point>> colorSquares[])
{
    squares.clear();

    Mat pyr, timg, gray0(image.size(), CV_8U), gray;

    // down-scale and upscale the image to filter out the noise
    pyrDown(image, pyr, Size(image.cols / 2, image.rows / 2));
    pyrUp(pyr, timg, image.size());
    vector<vector<Point>> contours;

    // find squares in every color plane of the image
    for (int c = 0; c < 3; ++c)
    {
        colorSquares[c].clear();

        int ch[] = {c, 0};
        mixChannels(&timg, 1, &gray0, 1, ch, 1);

        // try several threshold levels
        for (int l = 0; l < N; l++)
        {
            // hack: use Canny instead of zero threshold level.
            // Canny helps to catch squares with gradient shading
            if (l == 0)
            {
                // apply Canny. Take the upper threshold from slider
                // and set the lower to 0 (which forces edges merging)
                Canny(gray0, gray, 0, thresh, 5);
                // dilate canny output to remove potential
                // holes between edge segments
                dilate(gray, gray, Mat(), Point(-1, -1));
            }
            else
            {
                // apply threshold if l!=0:
                //     tgray(x,y) = gray(x,y) < (l+1)*255/N ? 255 : 0
                gray = gray0 >= (l + 1) * 255 / N;
            }

            // find contours and store them all as a list
            findContours(gray, contours, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);

            vector<Point> approx;

            // test each contour
            for (size_t i = 0; i < contours.size(); i++)
            {
                // approximate contour with accuracy proportional
                // to the contour perimeter
                approxPolyDP(Mat(contours[i]), approx, arcLength(Mat(contours[i]), true) * 0.02, true);

                // square contours should have 4 vertices after approximation
                // relatively large area (to filter out noisy contours)
                // and be convex.
                // Note: absolute value of an area is used because
                // area may be positive or negative - in accordance with the
                // contour orientation
                if (approx.size() == 4 &&
                    fabs(contourArea(Mat(approx))) > 1000 &&
                    isContourConvex(Mat(approx)))
                {
                    double maxCosine = 0;

                    for (int j = 2; j < 5; j++)
                    {
                        // find the maximum cosine of the angle between joint edges
                        double cosine = fabs(angle(approx[j % 4], approx[j - 2], approx[j - 1]));
                        maxCosine = MAX(maxCosine, cosine);
                    }

                    // if cosines of all angles are small
                    // (all angles are ~90 degree) then write quandrange
                    // vertices to resultant sequence
                    if (maxCosine < 0.3)
                        colorSquares[c].push_back(approx);
                }
            }
        }
    }

    FilterSquares(image, colorSquares, squares);
}

// the function draws all the squares in the image
static void DrawSquares(Mat &image, const vector<vector<Point>> &squares)
{
    for (size_t i = 0; i < squares.size(); i++)
    {
        const Point *p = &squares[i][0];
        int n = (int)squares[i].size();
        polylines(image, &p, &n, 1, true, Scalar(0, 255, 0), 3, CV_AA);
    }

    imshow(wndname, image);
}

// the function draws all the squares in the image on their respective color planes and colors
static void DrawSquares(Mat &image, const vector<vector<Point>> squares[])
{
    const int numOfColorPlanes = 3;
    Scalar color;
    const Scalar CV_BLUE(255, 0, 0);
    const Scalar CV_GREEN(0, 255, 0);
    const Scalar CV_RED(0, 0, 255);

    for (size_t c = 0; c < numOfColorPlanes; ++c)
    {
        switch (c)
        {
        case 0:
            color = CV_BLUE;
            break;
        case 1:
            color = CV_GREEN;
            break;
        case 2:
            color = CV_RED;
            break;
        default:
            break;
        }

        for (size_t i = 0; i < squares[c].size(); i++)
        {
            const Point *p = &squares[c][i][0];
            int n = (int)squares[c][i].size();
            polylines(image, &p, &n, 1, true, color, 3, CV_AA);
        }
    }

    imshow(wndname, image);
}

int main(int /*argc*/, char ** /*argv*/)
{
    vector<int> x(6);
    // Generate vector of ints from 1 to 6
    iota(begin(x), end(x), 1);

    vector<string> names(6);
    transform(x.cbegin(), x.cend(), names.begin(), [](const int in) { return "pic" + to_string(in) + ".png"; });

    help();
    namedWindow(wndname, 1);
    vector<vector<Point>> squares;
    vector<vector<Point>> colorSquares[3];

    for (auto name : names)
    {
        cout << "Name: " << name << endl;
        Mat image = imread(name, 1);
        if (image.empty())
        {
            cout << "Couldn't load " << name << endl;
            continue;
        }

        FindSquares(image, squares, colorSquares);
        DrawSquares(image, squares);
        //DrawSquares(image, colorSquares);

        int c = waitKey();
        if ((char)c == 27)
            break;
    }

    return 0;
}
