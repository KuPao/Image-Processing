#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <algorithm>

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene( void )
{
	
	ui_instance = Qt_Opengl_Framework::getInstance();
	
}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage( QString filePath )
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath )
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB( void )
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (! img_data )
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0 ; j < img_width ; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j*4), rgb + (out_offset + j*3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB( unsigned char *rgba, unsigned char *rgb )
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0 ; i < 3 ; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = i*img_width*3+j*3;
			int offset_rgba = i*img_width*4+j*4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k=0; k<3; k++)
				img_data[offset_rgba+k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			img_data[offset_rgba + rr] = (unsigned char)rgb[offset_rgb + rr] / 32 * 32;
			img_data[offset_rgba + gg] = (unsigned char)rgb[offset_rgb + gg] / 32 * 32;
			img_data[offset_rgba + bb] = (unsigned char)rgb[offset_rgb + bb] / 64 * 64;

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
struct GreaterCounts {
	typedef std::tuple<int, int, int> color;
	bool operator() (std::pair<color, int> const & a, std::pair<color, int> const & b) const {
		return a.second > b.second;
	}
};
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();

	typedef std::tuple<int, int, int> color;
	int factor = 8;
	std::map<color, int> colorCount;
	color currentColor;

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;

			rgb[offset_rgb + rr] = (rgb[offset_rgb + rr] / factor)*factor;
			rgb[offset_rgb + gg] = (rgb[offset_rgb + gg] / factor)*factor;
			rgb[offset_rgb + bb] = (rgb[offset_rgb + bb] / factor)*factor;
			
			currentColor = std::make_tuple(rgb[offset_rgb + rr], rgb[offset_rgb + gg], rgb[offset_rgb + bb]);

			if (colorCount.find(currentColor) == colorCount.end()) { //Not in map
				colorCount[color(currentColor)] = 1; //Initialize new color
			}
			else {
				colorCount[color(currentColor)]++;
			}
			
		}
	}

	if (colorCount.size() > 256) { //If there are fewer than 256 colors then finished

		std::vector<std::pair<color, int>> colors(colorCount.begin(), colorCount.end()); //Copy map key/value pairs to a vector
		sort(colors.begin(), colors.end(), GreaterCounts()); //Sort the vector by the value (color count) of the pair
		std::vector<color> reducedColors(256);

		int red, green, blue;
		for (int i = 0; i < 256; i++) {
			std::tie(red, green, blue) = colors[i].first;
			reducedColors[i] = std::make_tuple(red, green, blue);
		}

		int r0, g0, b0, r1, g1, b1;
		double minDis, dis;
		color closest;
		std::map <color, color> q_map;

		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgb = i * img_width * 3 + j * 3;
				int offset_rgba = i * img_width * 4 + j * 4;

				currentColor = std::make_tuple(rgb[offset_rgb + rr], rgb[offset_rgb + gg], rgb[offset_rgb + bb]);
				std::tie(r0, g0, b0) = currentColor;

				if (q_map.find(currentColor) == q_map.end()) { //Search for correct target color
				//Not in map, find closest of the 256
					minDis = INT_MAX; //Reset minimum distance for each source color
					for (int k = 0; k < 256; ++k) { //for each reduced color
						std::tie(r1, g1, b1) = reducedColors[k]; //unpack the tuple
						dis = sqrt(pow((double)(r1 - r0), 2) + pow((double)(g1 - g0), 2) + pow((double)(b1 - b0), 2)); //Calculate euclidean distance
						if (dis < minDis) {
							minDis = dis;
							closest = reducedColors[k];
						}
					}
					q_map[color(currentColor)] = closest; //Store closest color to the quantization map
				}

				std::tie(img_data[offset_rgba + rr], img_data[offset_rgba + gg], img_data[offset_rgba + bb]) = q_map[color(currentColor)]; 	//Overwrite the pixel color with the quantized color
			}
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			img_data[offset_rgba + rr] = img_data[offset_rgba + gg] = img_data[offset_rgba + bb] = gray > 127 ? WHITE : BLACK;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();

	srand(time(NULL));

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			float value/* = gray / 256 + (-0.2 + (((float)rand()) / (float)RAND_MAX) * 0.4)*/;
			value = (gray + (rand() % 128 - 64)) / 256;

			img_data[offset_rgba + rr] = img_data[offset_rgba + gg] = img_data[offset_rgba + bb] = value > 0.5 ? WHITE : BLACK;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	Gray();
	unsigned char *rgb = this->To_RGB();
	float* rgbf = new float[img_height * img_width * 3];
	//Gray();
	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset = i * img_width * 3 + j * 3;
			rgbf[offset + rr] = rgb[offset + rr];
			rgbf[offset + gg] = rgb[offset + gg];
			rgbf[offset + bb] = rgb[offset + bb];
		}
	}

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int right = i * img_width * 3 + (j + 1) * 3;
			int bottomLeft = (i + 1) * img_width * 3 + (j - 1) * 3;
			int bottom = (i + 1) * img_width * 3 + j * 3;
			int bottomRight = (i + 1) * img_width * 3 + (j + 1) * 3;

			float oldpixel = rgbf[offset_rgb + rr];
			float newpixel = oldpixel > 127 ? WHITE : BLACK;

			float error = oldpixel - newpixel;


			img_data[offset_rgba + rr] = img_data[offset_rgba + gg] = img_data[offset_rgba + bb] = (unsigned char)newpixel;
			img_data[offset_rgba + aa] = WHITE;

			if(j + 1 < img_width)
				rgbf[right + rr] = rgbf[right + gg] = rgbf[right + bb] += (float)7 / 16 * error;
			if (i + 1 < img_height) {
				if (j - 1 > 0)
					rgbf[bottomLeft + rr] = rgbf[bottomLeft + gg] = rgbf[bottomLeft + bb] += (float)3 / 16 * error;
				rgbf[bottom + rr] = rgbf[bottom + gg] = rgbf[bottom + bb] += (float)5 / 16 * error;
				if (j + 1 < img_width)
					rgbf[bottomRight + rr] = rgbf[bottomRight + gg] = rgbf[bottomRight + bb] += (float)1 / 16 * error;
			}
		}
		++i;
		if (i < img_height) {
			for (int j = img_width - 1; j >= 0; j--) {
				int offset_rgb = i * img_width * 3 + j * 3;
				int offset_rgba = i * img_width * 4 + j * 4;

				int right = i * img_width * 3 + (j - 1) * 3;
				int bottomLeft = (i + 1) * img_width * 3 + (j + 1) * 3;
				int bottom = (i + 1) * img_width * 3 + j * 3;
				int bottomRight = (i + 1) * img_width * 3 + (j - 1) * 3;

				float oldpixel = rgbf[offset_rgb + rr];
				float newpixel = oldpixel > 127 ? WHITE : BLACK;

				float error = oldpixel - newpixel;


				img_data[offset_rgba + rr] = img_data[offset_rgba + gg] = img_data[offset_rgba + bb] = (unsigned char)newpixel;
				img_data[offset_rgba + aa] = WHITE;

				if (j - 1 >= 0)
					rgbf[right + rr] = rgbf[right + gg] = rgbf[right + bb] += (float)7 / 16 * error;
				if (i + 1 < img_height) {
					if (j + 1 < img_width)
						rgbf[bottomLeft + rr] = rgbf[bottomLeft + gg] = rgbf[bottomLeft + bb] += (float)3 / 16 * error;
					rgbf[bottom + rr] = rgbf[bottom + gg] = rgbf[bottom + bb] += (float)5 / 16 * error;
					if (j - 1 >= 0)
						rgbf[bottomRight + rr] = rgbf[bottomRight + gg] = rgbf[bottomRight + bb] += (float)1 / 16 * error;
				}
			}
		}
	}

	delete[] rgb;
	delete[] rgbf;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
	std::vector<unsigned char> imgData;

	float sum = 0, avg = 0, ratio;
	int index, sort;

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			float gray = rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			sum += gray;
			imgData.push_back(gray);
		}
	}
	avg = sum / (img_width * img_height);
	std::sort(imgData.begin(), imgData.end());

	ratio = (double)avg / 255.0;
	index = imgData.size();
	sort = (1.0 - ratio)*index;

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			img_data[offset_rgba + rr] = img_data[offset_rgba + gg] = img_data[offset_rgba + bb] = gray > imgData[sort] ? WHITE : BLACK;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();

	float thresh_mat[4][4] = {{0.7059, 0.3529, 0.5882, 0.2353},
							  {0.0588, 0.9412, 0.8235, 0.4118},
							  {0.4706, 0.7647, 0.8824, 0.1176},
							  {0.1765, 0.5294, 0.2941, 0.6471}};

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			gray /= 256;

			float mask = thresh_mat[i % 4][j % 4];

			img_data[offset_rgba + rr] = img_data[offset_rgba + gg] = img_data[offset_rgba + bb] = gray > mask ? WHITE : BLACK;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
float Application::GetClosedRGB(float oldpixel, int rgb)
{

	return 0.0f;
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();
	float* rgbf = new float[img_height * img_width * 3];
	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset = i * img_width * 3 + j * 3;
			rgbf[offset + rr] = rgb[offset + rr];
			rgbf[offset + gg] = rgb[offset + gg];
			rgbf[offset + bb] = rgb[offset + bb];
		}
	}

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int rg[8] = { 0, 36, 72, 108, 144, 180, 216, 255 };
			int b[4] = { 0, 85, 170, 255 };
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int right = i * img_width * 3 + (j + 1) * 3;
			int bottomLeft = (i + 1) * img_width * 3 + (j - 1) * 3;
			int bottom = (i + 1) * img_width * 3 + j * 3;
			int bottomRight = (i + 1) * img_width * 3 + (j + 1) * 3;

			float oldR = rgbf[offset_rgb + rr], oldG = rgbf[offset_rgb + gg], oldB = rgbf[offset_rgb + bb];
			oldR = oldR > 255 ? 255 : (oldR < 0 ? 0 : oldR);
			oldG = oldG > 255 ? 255 : (oldG < 0 ? 0 : oldG);
			oldB = oldB > 255 ? 255 : (oldB < 0 ? 0 : oldB);
			float newR = rg[(int)oldR / 32], newG = rg[(int)oldG / 32], newB = b[(int)oldB / 64];

			img_data[offset_rgba + rr] = (unsigned char)newR;
			img_data[offset_rgba + gg] = (unsigned char)newG;
			img_data[offset_rgba + bb] = (unsigned char)newB;
			img_data[offset_rgba + aa] = WHITE;

			float errorR = oldR - newR, errorG = oldG - newG, errorB = oldB - newB;
			float error = sqrt(errorR * errorR + errorG * errorG + errorB * errorB);

			if (j + 1 < img_width) {
				rgbf[right + rr] += (float)7 / 16 * errorR;
				rgbf[right + gg] += (float)7 / 16 * errorG;
				rgbf[right + bb] += (float)7 / 16 * errorB;
				/*rgbf[right + rr] = rgbf[right + rr] > 255 ? 255 : (rgbf[right + rr] < 0 ? 0 : rgbf[right + rr]);
				rgbf[right + gg] = rgbf[right + gg] > 255 ? 255 : (rgbf[right + gg] < 0 ? 0 : rgbf[right + gg]);
				rgbf[right + bb] = rgbf[right + bb] > 255 ? 255 : (rgbf[right + bb] < 0 ? 0 : rgbf[right + bb]);*/
			}
			if (i + 1 < img_height) {
				if (j - 1 >= 0) {
					rgbf[bottomLeft + rr] += (float)3 / 16 * errorR;
					rgbf[bottomLeft + gg] += (float)3 / 16 * errorG;
					rgbf[bottomLeft + bb] += (float)3 / 16 * errorB;
					/*rgbf[bottomLeft + rr] = rgbf[bottomLeft + rr] > 255 ? 255 : (rgbf[bottomLeft + rr] < 0 ? 0 : rgbf[bottomLeft + rr]);
					rgbf[bottomLeft + gg] = rgbf[bottomLeft + gg] > 255 ? 255 : (rgbf[bottomLeft + gg] < 0 ? 0 : rgbf[bottomLeft + gg]);
					rgbf[bottomLeft + bb] = rgbf[bottomLeft + bb] > 255 ? 255 : (rgbf[bottomLeft + bb] < 0 ? 0 : rgbf[bottomLeft + bb]);*/
				}
				rgbf[bottom + rr] += (float)5 / 16 * errorR;
				rgbf[bottom + gg] += (float)5 / 16 * errorG;
				rgbf[bottom + bb] += (float)5 / 16 * errorB;
				/*rgbf[bottom + rr] = rgbf[bottom + rr] > 255 ? 255 : (rgbf[bottom + rr] < 0 ? 0 : rgbf[bottom + rr]);
				rgbf[bottom + gg] = rgbf[bottom + gg] > 255 ? 255 : (rgbf[bottom + gg] < 0 ? 0 : rgbf[bottom + gg]);
				rgbf[bottom + bb] = rgbf[bottom + bb] > 255 ? 255 : (rgbf[bottom + bb] < 0 ? 0 : rgbf[bottom + bb]);*/
				if (j + 1 < img_width) {
					rgbf[bottomRight + rr] += (float)1 / 16 * errorR; 
					rgbf[bottomRight + gg] += (float)1 / 16 * errorG;
					rgbf[bottomRight + bb] += (float)1 / 16 * errorB;
					/*rgbf[bottomRight + rr] = rgbf[bottomRight + rr] > 255 ? 255 : (rgbf[bottomRight + rr] < 0 ? 0 : rgbf[bottomRight + rr]);
					rgbf[bottomRight + gg] = rgbf[bottomRight + gg] > 255 ? 255 : (rgbf[bottomRight + gg] < 0 ? 0 : rgbf[bottomRight + gg]);
					rgbf[bottomRight + bb] = rgbf[bottomRight + bb] > 255 ? 255 : (rgbf[bottomRight + bb] < 0 ? 0 : rgbf[bottomRight + bb]);*/
				}
			}
		}
		if (++i < img_height) {
			for (int j = img_width - 1; j >= 0; j--) {
				int rg[8] = { 0, 36, 72, 108, 144, 180, 216, 255 };
				int b[4] = { 0, 85, 170, 255 };
				int offset_rgb = i * img_width * 3 + j * 3;
				int offset_rgba = i * img_width * 4 + j * 4;

				int right = i * img_width * 3 + (j - 1) * 3;
				int bottomLeft = (i + 1) * img_width * 3 + (j + 1) * 3;
				int bottom = (i + 1) * img_width * 3 + j * 3;
				int bottomRight = (i + 1) * img_width * 3 + (j - 1) * 3;

				float oldR = rgbf[offset_rgb + rr], oldG = rgbf[offset_rgb + gg], oldB = rgbf[offset_rgb + bb];
				oldR = oldR > 255 ? 255 : (oldR < 0 ? 0 : oldR);
				oldG = oldG > 255 ? 255 : (oldG < 0 ? 0 : oldG);
				oldB = oldB > 255 ? 255 : (oldB < 0 ? 0 : oldB);
				float newR = rg[(int)oldR / 32], newG = rg[(int)oldG / 32], newB = b[(int)oldB / 64];

				img_data[offset_rgba + rr] = (unsigned char)newR;
				img_data[offset_rgba + gg] = (unsigned char)newG;
				img_data[offset_rgba + bb] = (unsigned char)newB;
				img_data[offset_rgba + aa] = WHITE;

				float errorR = oldR - newR, errorG = oldG - newG, errorB = oldB - newB;

				if (j - 1 >= 0) {
					rgbf[right + rr] += (float)7 / 16 * errorR;
					rgbf[right + gg] += (float)7 / 16 * errorG;
					rgbf[right + bb] += (float)7 / 16 * errorB;
					/*rgbf[right + rr] = rgbf[right + rr] > 255 ? 255 : (rgbf[right + rr] < 0 ? 0 : rgbf[right + rr]);
					rgbf[right + gg] = rgbf[right + gg] > 255 ? 255 : (rgbf[right + gg] < 0 ? 0 : rgbf[right + gg]);
					rgbf[right + bb] = rgbf[right + bb] > 255 ? 255 : (rgbf[right + bb] < 0 ? 0 : rgbf[right + bb]);*/
				}
				if (i + 1 < img_height) {
					if (j + 1 < img_width) {
						rgbf[bottomLeft + rr] += (float)3 / 16 * errorR;
						rgbf[bottomLeft + gg] += (float)3 / 16 * errorG;
						rgbf[bottomLeft + bb] += (float)3 / 16 * errorB;
						/*rgbf[bottomLeft + rr] = rgbf[bottomLeft + rr] > 255 ? 255 : (rgbf[bottomLeft + rr] < 0 ? 0 : rgbf[bottomLeft + rr]);
						rgbf[bottomLeft + gg] = rgbf[bottomLeft + gg] > 255 ? 255 : (rgbf[bottomLeft + gg] < 0 ? 0 : rgbf[bottomLeft + gg]);
						rgbf[bottomLeft + bb] = rgbf[bottomLeft + bb] > 255 ? 255 : (rgbf[bottomLeft + bb] < 0 ? 0 : rgbf[bottomLeft + bb]);*/
					}
					rgbf[bottom + rr] += (float)5 / 16 * errorR;
					rgbf[bottom + gg] += (float)5 / 16 * errorG;
					rgbf[bottom + bb] += (float)5 / 16 * errorB;
					/*rgbf[bottom + rr] = rgbf[bottom + rr] > 255 ? 255 : (rgbf[bottom + rr] < 0 ? 0 : rgbf[bottom + rr]);
					rgbf[bottom + gg] = rgbf[bottom + gg] > 255 ? 255 : (rgbf[bottom + gg] < 0 ? 0 : rgbf[bottom + gg]);
					rgbf[bottom + bb] = rgbf[bottom + bb] > 255 ? 255 : (rgbf[bottom + bb] < 0 ? 0 : rgbf[bottom + bb]);*/
					if (j - 1 >= 0) {
						rgbf[bottomRight + rr] += (float)1 / 16 * errorR;
						rgbf[bottomRight + gg] += (float)1 / 16 * errorG;
						rgbf[bottomRight + bb] += (float)1 / 16 * errorB;
						/*rgbf[bottomRight + rr] = rgbf[bottomRight + rr] > 255 ? 255 : (rgbf[bottomRight + rr] < 0 ? 0 : rgbf[bottomRight + rr]);
						rgbf[bottomRight + gg] = rgbf[bottomRight + gg] > 255 ? 255 : (rgbf[bottomRight + gg] < 0 ? 0 : rgbf[bottomRight + gg]);
						rgbf[bottomRight + bb] = rgbf[bottomRight + bb] > 255 ? 255 : (rgbf[bottomRight + bb] < 0 ? 0 : rgbf[bottomRight + bb]);*/
					}
				}
			}
		}
	}

	delete[] rgb;
	delete[] rgbf;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering( double filter[][5] )
{
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::filtering( double **filter, int n )
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//apply filter
			for (int channel = 0; channel < 3; channel++) {
				double pixel = 0;

				for (int h = -n / 2; h < (float)n / 2; h++) {
					if (!(i + h > 0 && i + h < img_height))
						continue;
					for (int w = -n / 2; w < (float)n / 2; w++) {
						if (!(j + w > 0 && j + w < img_width))
							continue;

						int filter_rgb = (i + h) * img_width * 3 + (j + w) * 3;
						pixel += rgb[filter_rgb + channel] * filter[h + n / 2][w + n / 2];
					}
				}

				if (pixel > 255)
					pixel = 255;
				if (pixel < 0)
					pixel = 0;

				img_data[offset_rgba + channel] = pixel;
			}
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	double ** filter;
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			filter[i][j] = 0.04;

	filtering(filter, 5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	double ** filter;
	double bartlett[5][5] = {{1, 2, 3, 2, 1},
							 {2, 4, 6, 4, 2},
							 {3, 6, 9, 6, 3},
							 {2, 4, 6, 4, 2},
							 {1, 2, 3, 2, 1}};
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			filter[i][j] = bartlett[i][j] / 81;

	filtering(filter, 5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	double ** filter;
	double gaussian[5][5] = {{1,  4,  6,  4, 1},
							 {4, 16, 24, 16, 4},
							 {6, 24, 36, 24, 6},
							 {4, 16, 24, 16, 4},
							 {1,  4,  6,  4, 1}};
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			filter[i][j] = gaussian[i][j] / 256;

	filtering(filter, 5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N(unsigned int N )
{
	double ** filter;

	filter = (double **)malloc(N * sizeof(double *));
	for (int i = 0; i < N; i++)
		filter[i] = (double *)malloc(N * sizeof(double));

	int *binomial;
	binomial = (int *)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++) {
		// If this row not enough column for size, skip it.
		if (i + 1 < N)
			continue;
		// Calculate the coefficients.
		for (int j = 0, coef = 1; j <= i; j++) {
			coef = (!j || !i) ? 1 : (coef * (i - j + 1) / j);
			binomial[j] = coef;
		}
	}
	double sum = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == 0 || i == N - 1)
				filter[i][j] = binomial[j];
			else if (j == 0 || j == N - 1)
				filter[i][j] = binomial[i];
			else
				filter[i][j] = binomial[i] * binomial[j];
			sum += filter[i][j];
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			filter[i][j] = filter[i][j] / sum;

	filtering(filter, N);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	double ** filter;
	/*double edge[5][5] = {{-1,  -4,  -6,  -4, -1},
						 {-4, -16, -24, -16, -4},
						 {-6, -24, 220, -24, -6},
						 {-4, -16, -24, -16, -4},
						 {-1,  -4,  -6,  -4, -1}};*/
	double gaussian[5][5] = { {1,  4,  6,  4, 1},
							 {4, 16, 24, 16, 4},
							 {6, 24, 36, 24, 6},
							 {4, 16, 24, 16, 4},
							 {1,  4,  6,  4, 1} };

	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));

	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			filter[i][j] = gaussian[i][j] / -256;
	
	filter[2][2] += 1;

	filtering(filter, 5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();

	/*int N = 5;
	double ** filter;
	double edge[5][5];
	double enhance[5][5];

	edge[2][2] = enhance[2][2] = 1.0;

	filter = (double **)malloc(N * sizeof(double *));
	for (int i = 0; i < N; i++)
		filter[i] = (double *)malloc(N * sizeof(double));

	int *binomial;
	binomial = (int *)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++) {
		// If this row not enough column for size, skip it.
		if (i + 1 < N)
			continue;
		// Calculate the coefficients.
		for (int j = 0, coef = 1; j <= i; j++) {
			coef = (!j || !i) ? 1 : (coef * (i - j + 1) / j);
			binomial[j] = coef;
		}
	}
	double sum = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == 0 || i == N - 1)
				filter[i][j] = binomial[j];
			else if (j == 0 || j == N - 1)
				filter[i][j] = binomial[i];
			else
				filter[i][j] = binomial[i] * binomial[j];
			sum += filter[i][j];
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			filter[i][j] = filter[i][j] / sum;
			edge[i][j] = edge[i][j] - filter[i][j];
			enhance[i][j] += edge[i][j];
		}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			filter[i][j] = enhance[i][j];		

	filtering(filter, 5);*/

	double ** filter;
	double edge[5][5] = { {-1,  -4,  -6,  -4, -1},
						 {-4, -16, -24, -16, -4},
						 {-6, -24, 220, -24, -6},
						 {-4, -16, -24, -16, -4},
						 {-1,  -4,  -6,  -4, -1} };
	double enhance[5][5] = {{0, 0, 0, 0, 0},
							{0, 0, 0, 0, 0},
							{0, 0, 1, 0, 0},
							{0, 0, 0, 0, 0},
							{0, 0, 0, 0, 0}};
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			filter[i][j] = enhance[i][j] + edge[i][j] / 256;
	filtering(filter, 5);
	/*delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();*/
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	Resize(0.5);
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	Resize(2);
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

double Application::CubicInterpolation(double p0, double p1, double p2, double p3, double x)
{
	double a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
	double b = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3;
	double c = -0.5 * p0 + 0.5 * p2;
	double d = p1;
	return a * pow(x, 3) + b * pow(x, 2) + c * x + d;
}

double Application::BicubicInterpolation(int src[][4], double x, double y, int rgba)
{
	double pixel = 0;
	double coeffs[4][4];
	double v0, v1, v2, v3;

	coeffs[0][0] = img_data[src[1][1] + rgba];
	coeffs[0][1] = -.5*img_data[src[1][0] + rgba] + .5*img_data[src[1][2] + rgba];
	coeffs[0][2] = img_data[src[1][0] + rgba] - 2.5*img_data[src[1][1] + rgba] + 2 * img_data[src[1][2] + rgba] - .5*img_data[src[1][3] + rgba];
	coeffs[0][3] = -.5*img_data[src[1][0] + rgba] + 1.5*img_data[src[1][1] + rgba] - 1.5*img_data[src[1][2] + rgba] + .5*img_data[src[1][3] + rgba];
	coeffs[1][0] = -.5*img_data[src[0][1] + rgba] + .5*img_data[src[2][1] + rgba];
	coeffs[1][1] = .25*img_data[src[0][0] + rgba] - .25*img_data[src[0][2] + rgba] - .25*img_data[src[2][0] + rgba] + .25*img_data[src[2][2] + rgba];

	coeffs[1][2] = -.5*img_data[src[0][0] + rgba] + 1.25*img_data[src[0][1] + rgba] - img_data[src[0][2] + rgba] + .25*img_data[src[0][3] + rgba] +
		.5*img_data[src[2][0] + rgba] - 1.25*img_data[src[2][1] + rgba] + img_data[src[2][2] + rgba] - .25*img_data[src[2][3] + rgba];

	coeffs[1][3] = .25*img_data[src[0][0] + rgba] - .75*img_data[src[0][1] + rgba] + .75*img_data[src[0][2] + rgba] - .25*img_data[src[0][3] + rgba]
		- .25*img_data[src[2][0] + rgba] + .75*img_data[src[2][1] + rgba] - .75*img_data[src[2][2] + rgba] + .25*img_data[src[2][3] + rgba];

	coeffs[2][0] = img_data[src[0][1] + rgba] - 2.5*img_data[src[1][1] + rgba] + 2 * img_data[src[2][1] + rgba] - .5*img_data[src[3][1] + rgba];

	coeffs[2][1] = -.5*img_data[src[0][0] + rgba] + .5*img_data[src[0][2] + rgba] + 1.25*img_data[src[1][0] + rgba] - 1.25*img_data[src[1][2] + rgba]
		- img_data[src[2][0] + rgba] + img_data[src[2][2] + rgba] + .25*img_data[src[3][0] + rgba] - .25*img_data[src[3][2] + rgba];

	coeffs[2][2] = img_data[src[0][0] + rgba] - 2.5*img_data[src[0][1] + rgba] + 2 * img_data[src[0][2] + rgba] - .5*img_data[src[0][3] + rgba]
		- 2.5*img_data[src[1][0] + rgba] + 6.25*img_data[src[1][1] + rgba] - 5 * img_data[src[1][2] + rgba] + 1.25*img_data[src[1][3] + rgba] +
		2 * img_data[src[2][0] + rgba] - 5 * img_data[src[2][1] + rgba] + 4 * img_data[src[2][2] + rgba] - img_data[src[2][3] + rgba]
		- .5*img_data[src[3][0] + rgba] + 1.25*img_data[src[3][1] + rgba] - img_data[src[3][2] + rgba] + .25*img_data[src[3][3] + rgba];

	coeffs[2][3] = -.5*img_data[src[0][0] + rgba] + 1.5*img_data[src[0][1] + rgba] - 1.5*img_data[src[0][2] + rgba] + .5*img_data[src[0][3] + rgba] +
		1.25*img_data[src[1][0] + rgba] - 3.75*img_data[src[1][1] + rgba] + 3.75*img_data[src[1][2] + rgba] - 1.25*img_data[src[1][3] + rgba]
		- img_data[src[2][0] + rgba] + 3 * img_data[src[2][1] + rgba] - 3 * img_data[src[2][2] + rgba] + img_data[src[2][3] + rgba] +
		.25*img_data[src[3][0] + rgba] - .75*img_data[src[3][1] + rgba] + .75*img_data[src[3][2] + rgba] - .25*img_data[src[3][3] + rgba];

	coeffs[3][0] = -.5*img_data[src[0][1] + rgba] + 1.5*img_data[src[1][1] + rgba] - 1.5*img_data[src[2][1] + rgba] + .5*img_data[src[3][1] + rgba];
	coeffs[3][1] = .25*img_data[src[0][0] + rgba] - .25*img_data[src[0][2] + rgba] 
		- .75*img_data[src[1][0] + rgba] + .75*img_data[src[1][2] + rgba] +
		.75*img_data[src[2][0] + rgba] - .75*img_data[src[2][2] + rgba] 
		- .25*img_data[src[3][0] + rgba] + .25*img_data[src[3][2] + rgba];
	coeffs[3][2] = -.5*img_data[src[0][0] + rgba] + 1.25*img_data[src[0][1] + rgba] - img_data[src[0][2] + rgba] + .25*img_data[src[0][3] + rgba] +
		1.5*img_data[src[1][0] + rgba] - 3.75*img_data[src[1][1] + rgba] + 3 * img_data[src[1][2] + rgba] - .75*img_data[src[1][3] + rgba]
		- 1.5*img_data[src[2][0] + rgba] + 3.75*img_data[src[2][1] + rgba] - 3 * img_data[src[2][2] + rgba] + .75*img_data[src[2][3] + rgba] +
		.5*img_data[src[3][0] + rgba] - 1.25*img_data[src[3][1] + rgba] + img_data[src[3][2] + rgba] - .25*img_data[src[3][3] + rgba];
	coeffs[3][3] = .25*img_data[src[0][0] + rgba] - .75*img_data[src[0][1] + rgba] + .75*img_data[src[0][2] + rgba] - .25*img_data[src[0][3] + rgba]
		- .75*img_data[src[1][0] + rgba] + 2.25*img_data[src[1][1] + rgba] - 2.25*img_data[src[1][2] + rgba] + .75*img_data[src[1][3] + rgba] +
		.75*img_data[src[2][0] + rgba] - 2.25*img_data[src[2][1] + rgba] + 2.25*img_data[src[2][2] + rgba] - .75*img_data[src[2][3] + rgba]
		- .25*img_data[src[3][0] + rgba] + .75*img_data[src[3][1] + rgba] - .75*img_data[src[3][2]] + rgba + .25*img_data[src[3][3] + rgba];

	double x2 = x * x;
	double x3 = x2 * x;
	double y2 = y * y;
	double y3 = y2 * y;

	pixel = (coeffs[0][0] + coeffs[0][1] * y + coeffs[0][2] * y2 + coeffs[0][3] * y3) +
		(coeffs[1][0] + coeffs[1][1] * y + coeffs[1][2] * y2 + coeffs[1][3] * y3) * x +
		(coeffs[2][0] + coeffs[2][1] * y + coeffs[2][2] * y2 + coeffs[2][3] * y3) * x2 +
		(coeffs[3][0] + coeffs[3][1] * y + coeffs[3][2] * y2 + coeffs[3][3] * y3) * x3;

	/*v0 = CubicInterpolation(img_data[src[0][0] + rgba], img_data[src[0][1] + rgba], img_data[src[0][2] + rgba], img_data[src[0][3] + rgba], x);
	v1 = CubicInterpolation(img_data[src[1][0] + rgba], img_data[src[1][1] + rgba], img_data[src[1][2] + rgba], img_data[src[1][3] + rgba], x);
	v2 = CubicInterpolation(img_data[src[2][0] + rgba], img_data[src[2][1] + rgba], img_data[src[2][2] + rgba], img_data[src[2][3] + rgba], x);
	v3 = CubicInterpolation(img_data[src[3][0] + rgba], img_data[src[3][1] + rgba], img_data[src[3][2] + rgba], img_data[src[3][3] + rgba], x);
	pixel = CubicInterpolation(v0, v1, v2, v3, y);
	if (pixel != 0)
		pixel = pixel;*/
	return pixel > 255 ? 255 : (pixel < 0 ? 0 : pixel);
}

double Application::BilinearInterpolation(int src[][4], int srcX[4], int srcY[4], double x, double y, int rgba)
{
	double factor1;
	double factor2;
	double factor3;
	double factor4;
	if (srcX[2] == srcX[1]) // avoid divide by zero
	{
		factor1 = 1; // force calculatione to one point
		factor2 = 0;
	}
	else
	{
		factor1 = (((double)srcX[2] - x) / ((double)srcX[2] - (double)srcX[1]));
		factor2 = ((x - (double)srcX[1]) / ((double)srcX[2] - (double)srcX[1]));
	}

	double R1 = factor1 * (double)img_data[src[1][1] + rgba] + factor2 * (double)img_data[src[1][2] + rgba];
	double R2 = factor1 * (double)img_data[src[2][1] + rgba] + factor2 * (double)img_data[src[2][2] + rgba];

	if (srcY[2] == srcY[1])
	{
		factor3 = 1;
		factor4 = 0;
	}
	else
	{
		factor3 = ((double)srcY[2] - y) / ((double)srcY[2] - (double)srcY[1]);
		factor4 = (y - (double)srcY[1]) / ((double)srcY[2] - (double)srcY[1]);
	}
	double pixel = (unsigned int)((factor3 * R1) + (factor4*R2));
	
	return pixel > 255 ? 255 : (pixel < 0 ? 0 : pixel);
}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize( float scale )
{
	unsigned char *rgb = this->To_RGB();

	int new_img_width = scale * img_width, new_img_height = scale * img_height;
	unsigned char* dest_img = new unsigned char[new_img_width * new_img_height * 4];
	//img_data = new unsigned char[new_img_width * new_img_height * 4];

	for (int i = 0; i < new_img_height; i++) {
		for (int j = 0; j < new_img_width; j++) {
			int offset_rgba = i * new_img_width * 4 + j * 4;
			double x = (double)j / scale, y = (double)i / scale;
			
			int srcX[4], srcY[4];
			srcX[1] = (int)floor(x);
			srcX[0] = std::max(0, srcX[1] - 1);
			srcX[2] = std::min(img_width - 1, srcX[1] + 1);
			srcX[3] = std::min(img_width - 1, srcX[1] + 2);
			srcY[1] = (int)floor(y);
			srcY[0] = std::max(0, srcY[1] - 1);
			srcY[2] = std::min(img_height - 1, srcY[1] + 1);
			srcY[3] = std::min(img_height - 1, srcY[1] + 2);
			
			int src[4][4];
			for (int a = 0; a < 4; a++) 
				for (int b = 0; b < 4; b++)
					src[a][b] = srcY[a] * img_width * 4 + srcX[b] * 4;

			for (int rgba = 0; rgba < 4; rgba++)
				dest_img[offset_rgba + rgba] = (unsigned char)BilinearInterpolation(src, srcX, srcY, x, y, rgba);
				//dest_img[offset_rgba + rgba] = (unsigned char)img_data[src[1][1] + rgba];
		}
	}

	img_data = dest_img;
	img_width = new_img_width;
	img_height = new_img_height;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate( float angleDegrees )
{
	int times = (int)angleDegrees / 360;
	angleDegrees = angleDegrees - 360 * times;
	float radian = -(angleDegrees * M_PI / 180);

	double cosA = cos(radian);
	double sinA = sin(radian);
	int w = img_width;
	int h = img_height;
	int newW = round(w * fabs(cosA) + h * fabs(sinA));
	int newH = round(h * fabs(cosA) + w	 * fabs(sinA));
	newW = w;
	newH = h;

	double affine[] = { cosA, sinA, (1 - cosA) * w / 2. - sinA * h / 2.,
					   -sinA, cosA, sinA * w / 2. + (1 - cosA) * h / 2. };

	unsigned char* dest_img = new unsigned char[newW * newH * 4];
	for (int y = 0; y < newH; ++y)
	{
		for (int x = 0; x < newW; ++x)
		{
			int offset_rgba = y * newW * 4 + x * 4;

			int X = x * affine[0] + y * affine[1] + affine[2];
			int Y = x * affine[3] + y * affine[4] + affine[5];

			/*double X = x * affine[0] + y * affine[1] + affine[2], Y = x * affine[3] + y * affine[4] + affine[5];

			int srcX[4], srcY[4];
			srcX[1] = (int)floor(x);
			srcX[0] = std::max(0, srcX[1] - 1);
			srcX[2] = std::min(img_width - 1, srcX[1] + 1);
			srcX[3] = std::min(img_width - 1, srcX[1] + 2);
			srcY[1] = (int)floor(y);
			srcY[0] = std::max(0, srcY[1] - 1);
			srcY[2] = std::min(img_height - 1, srcY[1] + 1);
			srcY[3] = std::min(img_height - 1, srcY[1] + 2);

			int src[4][4];
			for (int a = 0; a < 4; a++)
				for (int b = 0; b < 4; b++)
					src[a][b] = srcY[a] * img_width * 4 + srcX[b] * 4;*/

			if ((unsigned)X < w && (unsigned)Y < h)
			{
				int old_offset = Y * w * 4 + X * 4;
				for (int rgba = 0; rgba < 4; rgba++)
					dest_img[offset_rgba + rgba] = img_data[old_offset + rgba];
				/*for (int rgba = 0; rgba < 4; rgba++)
					dest_img[offset_rgba + rgba] = (unsigned char)BilinearInterpolation(src, srcX, srcY, x, y, rgba);*/
			}
			else {
				dest_img[offset_rgba + rr] = WHITE, dest_img[offset_rgba + gg] = WHITE, dest_img[offset_rgba + bb] = WHITE, dest_img[offset_rgba + aa] = BLACK;
			}
		}
	}

	img_data = dest_img;
	img_width = newW;
	img_height = newH;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge( QString filePath )
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image( int tMethod )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	// new canvas
	unsigned char* dst_img = this->To_RGB();

	int brushes[3];
	brushes[0] = 7;
	brushes[1] = 3;
	brushes[2] = 1;
	for (int i = 0; i < 3; i++) {
		unsigned char* reference = new unsigned char[img_height * img_width * 3];
		memcpy(reference, dst_img, img_height * img_width * 3);
		NPR_Gaussian_N(2 * brushes[i] + 1, reference);

		NPR_Paint_Layer(img_data, reference, brushes[i]);
	}
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::NPR_Paint_Layer( unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize )
{
	std::vector<Stroke> strokes;
	float* diff = new float[img_width * img_height];
	difference(tCanvas, tReferenceImage, diff);

	int grid = static_cast<int>(ceil(gridSizeFactor * tBrushSize));

	srand((unsigned)time(NULL));

	for (int y = 0; y < img_height; y += grid)
		for (int x = 0; x < img_width; x += grid) {

			float areaError = 0;

			for(int i = floor((float)grid * -1 / 2); i < ceil((float)grid / 2); i++)
				for (int j = floor((float)grid * -1 / 2); j < ceil((float)grid / 2); j++) {
					int v = y + i;
					int u = x + j;
					if (v < 0 || v >= img_height || u < 0 || u >= img_width)
						continue;
					else
						areaError += diff[v * img_width + u];
				}

			areaError / grid / grid;

			if (areaError > 25) {
				float max = 0;
				int x1, y1;
				for (int i = floor((float)grid * -1 / 2); i < ceil((float)grid / 2); i++)
					for (int j = floor((float)grid * -1 / 2); j < ceil((float)grid / 2); j++) {
						int v = y + i;
						int u = x + j;
						if (v < 0 || v >= img_height || u < 0 || u >= img_width)
							continue;
						else if (max < diff[v * img_width + u]) {
							max = diff[v * img_width + u];
							x1 = u;
							y1 = v;
						}
							
					}

				int r, g, b;
				int offset_rgb = y1 * img_width * 3 + x1 * 3;
				r = tReferenceImage[offset_rgb + rr];
				g = tReferenceImage[offset_rgb + gg];
				b = tReferenceImage[offset_rgb + bb];

				strokes.push_back(Stroke(tBrushSize / 2, x1, y1, r, g, b, WHITE));
			}
		}

	for (; strokes.size() > 0;) {
		int id = rand() % strokes.size();
		Paint_Stroke(strokes[id]);
		strokes.erase(strokes.begin() + id);
	}
}

void Application::difference(unsigned char *img1, unsigned char *img2, float* diff)
{
	for (unsigned int row = 0; row < img_height; ++row)
		for (unsigned int col = 0; col < img_width; ++col) {
			int offset_rgb = row * img_width * 3 + col * 3;
			int offset_rgba = row * img_width * 4 + col * 4;

			float r, g, b;
			r = img1[offset_rgba + rr] - img2[offset_rgb + rr];
			g = img1[offset_rgba + gg] - img2[offset_rgb + gg];
			b = img1[offset_rgba + bb] - img2[offset_rgb + bb];

			diff[row * img_width + col] = sqrt(r*r + g*g + b*b);
		}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke( const Stroke& s )
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) 
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) 
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height)) 
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared) 
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				} 
				else if (dist_squared == radius_squared + 1) 
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}

void Application::NPR_filtering(double **filter, int n, unsigned char* dst)
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//apply filter
			for (int channel = 0; channel < 3; channel++) {
				double pixel = 0;

				for (int h = -n / 2; h < (float)n / 2; h++) {
					if (!(i + h > 0 && i + h < img_height))
						continue;
					for (int w = -n / 2; w < (float)n / 2; w++) {
						if (!(j + w > 0 && j + w < img_width))
							continue;

						int filter_rgb = (i + h) * img_width * 3 + (j + w) * 3;
						pixel += rgb[filter_rgb + channel] * filter[h + n / 2][w + n / 2];
					}
				}

				if (pixel > 255)
					pixel = 255;
				if (pixel < 0)
					pixel = 0;

				dst[offset_rgb + channel] = pixel;
			}
		}
	}
	delete[] rgb;
}

void Application::NPR_Gaussian_N(unsigned int N, unsigned char* dst)
{
	double ** filter;

	filter = (double **)malloc(N * sizeof(double *));
	for (int i = 0; i < N; i++)
		filter[i] = (double *)malloc(N * sizeof(double));

	int *binomial;
	binomial = (int *)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++) {
		// If this row not enough column for size, skip it.
		if (i + 1 < N)
			continue;
		// Calculate the coefficients.
		for (int j = 0, coef = 1; j <= i; j++) {
			coef = (!j || !i) ? 1 : (coef * (i - j + 1) / j);
			binomial[j] = coef;
		}
	}
	double sum = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == 0 || i == N - 1)
				filter[i][j] = binomial[j];
			else if (j == 0 || j == N - 1)
				filter[i][j] = binomial[i];
			else
				filter[i][j] = binomial[i] * binomial[j];
			sum += filter[i][j];
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			filter[i][j] = filter[i][j] / sum;

	NPR_filtering(filter, N, dst);
}



///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}





