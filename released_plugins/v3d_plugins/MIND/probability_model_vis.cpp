/**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */
#include "probability_model_vis.h"

/**
 * @brief Function to map a value between 0 and 1 to a color in the colormap
 *
 * @param value - the value to map (between 0 and 1)
 */
std::tuple<unsigned char, unsigned char, unsigned char> colormap(double value) {
  // Ensure value is between 0 and 1
  value = std::clamp(value, 0.0, 1.0);
  int index = static_cast<int>(value * (INFERNO_COLORMAP_SIZE - 1));

  double r = std::clamp(std::get<0>(INFERNO_COLORMAP[index]), 0.0, 1.0);
  double g = std::clamp(std::get<1>(INFERNO_COLORMAP[index]), 0.0, 1.0);
  double b = std::clamp(std::get<2>(INFERNO_COLORMAP[index]), 0.0, 1.0);

  return {static_cast<unsigned char>(r * 255),
          static_cast<unsigned char>(g * 255),
          static_cast<unsigned char>(b * 255)};
}

/**
 * @brief Function to visualize PCA results
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 */
void visualizeProbabilityModel_func(V3DPluginCallback2 &callback,
                                    QWidget *parent) {
  // Load data
  QString filename = QFileDialog::getOpenFileName(
      parent, "Open Probabilistic Model", "", "Binary Files (*.bin)");

  if (filename.isEmpty()) {
    printf("No file selected.\n");
    return;
  }

  std::vector<double> data;
  V3DLONG dim_X, dim_Y, dim_Z;
  cellSegmentation::class_segmentationMain::loadProbabilityModel(
      filename.toStdString().c_str(), data, dim_X, dim_Y, dim_Z);

  Image4DSimple *p4DImage = new Image4DSimple();
  p4DImage->createBlankImage(dim_X, dim_Y, dim_Z, 3, V3D_UINT8);

  double min = data[0];
  double max = data[0];

  for (V3DLONG i = 0; i < dim_X * dim_Y * dim_Z; i++) {
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  unsigned char *pixels = p4DImage->getRawData();
  int channelSize = dim_X * dim_Y * dim_Z;
  for (V3DLONG i = 0; i < dim_X * dim_Y * dim_Z; i++) {
    auto [r, g, b] = colormap((data[i] - min) / (max - min));
    pixels[i] = r;
    pixels[i + channelSize] = g;
    pixels[i + 2 * channelSize] = b;
  }

  v3dhandle newwin = callback.newImageWindow("Probabilistic Model");
  callback.setImage(newwin, p4DImage);

  // open legend
  probabilityModelLegend_func(callback, parent);
}

/**
 * @brief Function to create a legend for the probabilistic model
 */
void probabilityModelLegend_func(V3DPluginCallback2 &callback,
                                 QWidget *parent) {
  const int inner_width = 512;  // Width of the color bar
  const int inner_height = 60;  // Height of the color bar
  const int padding_x = 40;     // Padding around the color bar
  const int padding_y = 15;
  const int full_width = inner_width + 2 * padding_x;
  const int full_height = inner_height + padding_y;

  QImage legendImage(full_width, full_height, QImage::Format_ARGB32);
  legendImage.fill(Qt::white);

  // Draw the colormap in the top portion (say the top 20 pixels)
  const int colorBarHeight = 30;
  for (int x = padding_x; x < inner_width + padding_x; ++x) {
    double value = double(x - padding_x) / (inner_width - 1);  // [0,1]
    auto [r, g, b] = colormap(value);
    QColor color(static_cast<int>(r), static_cast<int>(g), static_cast<int>(b));

    // Fill a small vertical column of color
    for (int y = 0; y < colorBarHeight; ++y) {
      legendImage.setPixelColor(x, y, color);
    }
  }

  // Create a QPainter to draw numeric ticks/labels
  QPainter painter(&legendImage);
  painter.setPen(Qt::black);
  painter.setFont(QFont("Arial", 12));

  // Define which ticks to draw. Here we do 0.0, 0.25, 0.5, 0.75, 1.0
  QList<double> ticks = {0.0, 0.25, 0.5, 0.75, 1.0};
  for (double t : ticks) {
    int xPos = int(t * (inner_width - 1) + padding_x);
    // Vertical position below the color bar
    int textY = colorBarHeight + 32;
    // Draw the numeric label
    QString label = QString::number(t, 'f', 2);  // e.g. "0.00", "0.25", ...
    painter.drawText(xPos - 25, textY, label);
    // Optionally, draw a small tick line at each label
    painter.drawLine(xPos, colorBarHeight, xPos, colorBarHeight + 10);
  }

  // Set up a dialog to display our legend
  QDialog *subWindow = new QDialog(callback.getVaa3DMainWindow());
  subWindow->setWindowTitle("Colormap Legend");
  subWindow->setWindowFlags(
      Qt::Window | Qt::WindowTitleHint | Qt::CustomizeWindowHint |
      Qt::WindowCloseButtonHint | Qt::MSWindowsFixedSizeDialogHint);

  QLabel *label = new QLabel(subWindow);
  label->setPixmap(QPixmap::fromImage(legendImage));
  label->setAlignment(Qt::AlignCenter);
  // If you do not want the label to stretch, omit "setScaledContents(true)"
  label->setScaledContents(true);

  QVBoxLayout *layout = new QVBoxLayout(subWindow);
  layout->addWidget(label);
  subWindow->setLayout(layout);

  // Adjust as needed for an initial display size
  subWindow->resize(600, 120);
  subWindow->show();
}