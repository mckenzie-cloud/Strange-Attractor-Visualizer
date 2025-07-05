
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>

#include <SFML/Graphics.hpp>
#include "SOLVER/Solver.h"
#include "Attractor.h"

#define VISUALIZER_WINDOW_SIZE 640
#define PLOTTER_WINDOW_SIZE 480

#define N_DIMENSION 1    // Don't modify. Use for rotation and projection matrix
#define M_DIMENSION 3    // Don't modify. Use for rotation and projection matrix
#define ODE_DIMENSION 3  // Don't modify. We are solving for 3-dimensional ODEs.

struct Quaternion {
    long double w, x, y, z;

    // w, x, y, z
    Quaternion(long double w, long double x, long double y, long double z) : w(w), x(x), y(y), z(z) 
    {}

    Quaternion operator*(const Quaternion &q) const {

        long double ww, xx, yy, zz;
        ww = w*q.w - x*q.x - y*q.y - z*q.z;
        xx = w*q.x + x*q.w + y*q.z - z*q.y;
        yy = w*q.y - x*q.z + y*q.w + z*q.x;
        zz = w*q.z + x*q.y - y*q.x + z*q.w;

        return Quaternion(ww, xx, yy, zz);
    }

    void conjugate() {       // Quaternion Conjugate (w - ix - jy - kz)
        x = -x;
        y = -y;
        z = -z;
    }
};

void visualizeODE(std::vector<std::array<long double, 2>> &projectedPoints,
            sf::RenderWindow &window, float visualizingScalingFactor);

void plotFunction(std::vector<std::vector<long double>> &y_out,
            std::vector<long double> &t_out, std::vector<long double> &maxPoints,
            sf::RenderWindow &targetWindow, float offsetX, float visualizingScalingFactor, 
            float lineHeight, int targetAxis);

long double vectorMag3(std::array<long double, ODE_DIMENSION> &vec);
void normalize3D(std::array<long double, ODE_DIMENSION> &vec);
void qRotation(Quaternion &p, Quaternion &q, Quaternion &p_prime);
Quaternion axisAngleToQuaternion(std::array<long double, ODE_DIMENSION> &axis, long double angleRadians);

std::array<std::string, 10> attractors = {       
    "LORENZ",
    "THOMAS",
    "AIZAWA",
    "DADRAS",
    "ROSSLER",
    "SPROTTB",
    "ARNEODO",
    "LORENZ83",
    "HALVORSEN",
    "SPROTT_LINZF"
};

std::array<std::string, ODE_DIMENSION> axis = {"X-axis", "Y-axis", "Z-axis"};

// Color
sf::Color axisLineColor(237, 53, 0);
sf::Color plottingWindowBGColor(255, 252, 251);
sf::Color plottingLineColor(255, 180, 180);
sf::Color plottingWindowTextColor(138, 0, 0);

int main(void)
{
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow visualizationWindow(sf::VideoMode(VISUALIZER_WINDOW_SIZE, VISUALIZER_WINDOW_SIZE), "Visualizing Attractor | Copr. MacKenzie 2025", sf::Style::Close, settings);
    sf::RenderWindow plottingWindow(sf::VideoMode(PLOTTER_WINDOW_SIZE, PLOTTER_WINDOW_SIZE), "Plotting Function", sf::Style::Close, settings);
    visualizationWindow.setVerticalSyncEnabled(true);   // Activating vertical synchronization will limit the number of frames displayed to the refresh rate of the monitor.

    long double t0 = 0.0;
    long double tFinal;
    std::vector<long double> y0(ODE_DIMENSION);
    std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> targetFunction;

    Solver *solver = nullptr;
    long double absTol = 1e-6, relTol = 1e-6;
    StepSizeController::Controllers controller = StepSizeController::STANDARD; // Choose a controller

    std::vector<long double> limitPoints(ODE_DIMENSION);

    float visualizingScalingFactor = 5.0f;
    float lineHeight = 100.f;
    float plottingWindowSize = static_cast<float>(PLOTTER_WINDOW_SIZE - 80);
    float plottingWindowOffsetX = static_cast<float>((PLOTTER_WINDOW_SIZE - plottingWindowSize) / 2);
    float plottingScalingFactor = 0.f;

    sf::RectangleShape axis1, axis2, maxPointLine, originPointLine, tFinalLine;

    axis1.setPosition(sf::Vector2f(plottingWindowOffsetX, static_cast<float>(PLOTTER_WINDOW_SIZE/2)));
    axis1.setSize(sf::Vector2f(plottingWindowSize, 3.f));
    axis1.setFillColor(axisLineColor);

    axis2.setPosition(sf::Vector2f(plottingWindowOffsetX - 3.f, static_cast<float>(PLOTTER_WINDOW_SIZE/2)));
    axis2.setSize(sf::Vector2f(2.f * lineHeight, 3.f));
    axis2.setFillColor(axisLineColor);
    axis2.rotate(-90.f);

    maxPointLine.setPosition(sf::Vector2f(plottingWindowOffsetX - 6.f, static_cast<float>(PLOTTER_WINDOW_SIZE/2 - 1.5f) - 2.f * lineHeight));
    maxPointLine.setSize(sf::Vector2f(10.f, 1.5f));
    maxPointLine.setFillColor(axisLineColor);

    originPointLine.setPosition(sf::Vector2f(plottingWindowOffsetX - 6.f, static_cast<float>(PLOTTER_WINDOW_SIZE/2 - 1.5f) - lineHeight));
    originPointLine.setSize(sf::Vector2f(10.f, 1.5f));
    originPointLine.setFillColor(axisLineColor);

    tFinalLine.setPosition(sf::Vector2f(plottingWindowOffsetX + plottingWindowSize, static_cast<float>(PLOTTER_WINDOW_SIZE/2 + 6.f)));
    tFinalLine.setSize(sf::Vector2f(10.f, 1.5f));
    tFinalLine.setFillColor(axisLineColor);
    tFinalLine.rotate(-90.f);

    sf::Vector2f attractorListBtnDim(80.f, 16.f);
    sf::Vector2f axisListBtnDim(80.f, 16.f);
    sf::Vector2f attractorListBtnPos(VISUALIZER_WINDOW_SIZE - (attractorListBtnDim.x + 5.f), 5.f);
    sf::Vector2f axisListBtnPos(PLOTTER_WINDOW_SIZE/2 - 40.f, (PLOTTER_WINDOW_SIZE/2 + 8.f) + 100.f);

    // Load texture
    sf::Texture texture;
    texture.loadFromFile("../assets/icon/btn_texture.png");

    // Load sprite
    sf::Sprite attractor_list_btn_sprite(texture);
    attractor_list_btn_sprite.setPosition(attractorListBtnPos);

    sf::Sprite axis_list_btn_sprite(texture);
    axis_list_btn_sprite.setPosition(axisListBtnPos);

    // Load font
    sf::Font font;
    font.loadFromFile("../assets/fonts/PTSans-Regular.ttf");

    // Initialize Texts
    sf::Text messageText1("", font, 20);
    sf::Text messageText2("Use [Up/Down] arrow key to zoom [In/Out].", font, 20);
    sf::Text axisListText("", font, 12);
    sf::Text axis1LabelText("", font, 20);
    sf::Text axis2LabelText("Time", font, 20);
    sf::Text tFinalLabelText("", font, 12);
    sf::Text minPointLabelText("", font, 12);
    sf::Text maxPointLabelText("", font, 12);
    sf::Text attractorListText("", font, 12);
    sf::Text originPointLabelText("0", font, 12);
    sf::Text axisListbtnLabelText("Change axis", font, 12);
    sf::Text attractorListBtnLabelText("Select an Attractor: ", font, 12);


    // Set text color
    axisListText.setFillColor(plottingWindowTextColor);
    axis1LabelText.setFillColor(plottingWindowTextColor);
    axis2LabelText.setFillColor(plottingWindowTextColor);
    tFinalLabelText.setFillColor(plottingWindowTextColor);
    minPointLabelText.setFillColor(plottingWindowTextColor);
    maxPointLabelText.setFillColor(plottingWindowTextColor);
    originPointLabelText.setFillColor(plottingWindowTextColor);
    axisListbtnLabelText.setFillColor(plottingWindowTextColor);

    // Set Text positions
    float messageText2Width = messageText2.getLocalBounds().width;
    float messageText2Height = messageText2.getLocalBounds().height;
    float axis2LabelTextWidth = axis2LabelText.getLocalBounds().width;
    float originPointLabelTextWidth = originPointLabelText.getLocalBounds().width;
    float axisListbtnLabelTextWidth = axisListbtnLabelText.getLocalBounds().width;
    float attractorListBtnLabelTextWidth = attractorListBtnLabelText.getLocalBounds().width;
    float originPointLabelTextHeight = originPointLabelText.getLocalBounds().height;
    float axisListbtnLabelTextHeight = axisListbtnLabelText.getLocalBounds().height;
    float attractorListBtnLabelTextHeight = attractorListBtnLabelText.getLocalBounds().height;

    float messageText2offSetX = (VISUALIZER_WINDOW_SIZE - messageText2Width) / 2.f;
    float axisListBtnLabelTextOffsetX = (axisListBtnDim.x - axisListbtnLabelTextWidth) / 2.f;

    sf::Vector2f attractorListTextPos, axisListTextPos;
    messageText2.setPosition(sf::Vector2f(messageText2offSetX, VISUALIZER_WINDOW_SIZE - messageText2Height - 10.f));
    axis2LabelText.setPosition(sf::Vector2f(plottingWindowOffsetX + (plottingWindowSize/2.f - axis2LabelTextWidth/2.f), PLOTTER_WINDOW_SIZE/2.f));
    tFinalLabelText.setPosition(sf::Vector2f(plottingWindowOffsetX + plottingWindowSize, (PLOTTER_WINDOW_SIZE/2.f) + 10.f));
    minPointLabelText.setPosition(sf::Vector2f(0.f, (PLOTTER_WINDOW_SIZE/2.f)));
    originPointLabelText.setPosition(sf::Vector2f(plottingWindowOffsetX - 2.f * originPointLabelTextWidth, (PLOTTER_WINDOW_SIZE/2.f - (lineHeight + originPointLabelTextHeight))));
    axisListbtnLabelText.setPosition(sf::Vector2f(axisListBtnPos.x + axisListBtnLabelTextOffsetX, axisListBtnPos.y - axisListbtnLabelTextHeight - 5.f));
    attractorListBtnLabelText.setPosition(sf::Vector2f(attractorListBtnPos.x - attractorListBtnLabelTextWidth, attractorListBtnLabelTextHeight/2.f));


    // Rotate text
    axis1LabelText.rotate(-90.f);

    // default value for our state machine
    bool state_visualize = false, state_solveSystem = false;
    bool state_isAttractorListBtnClicked = false, state_isAxisListBtnClicked = false;

    int attractorId = -1;
    int targetAxis  = 0;        // default is x-axis

    std::array<long double, ODE_DIMENSION> rotationAxis = {0.70710678118, 0.0, 0.70710678118};
    long double theta = 0.0;

    while (visualizationWindow.isOpen() && plottingWindow.isOpen())
    {
        // clear background
        visualizationWindow.clear(sf::Color::Black);
        plottingWindow.clear(plottingWindowBGColor);
        
        sf::Event visualizationWindowEvent;
        while (visualizationWindow.pollEvent(visualizationWindowEvent))
        {
            if (visualizationWindowEvent.type == sf::Event::Closed)
            {
                visualizationWindow.close();
            }

            // handle keyboard events
            if (visualizationWindowEvent.type == sf::Event::KeyPressed) {
                
                if (visualizationWindowEvent.key.code == sf::Keyboard::Up) {
                    visualizingScalingFactor += 2;
                }
                if (visualizationWindowEvent.key.code == sf::Keyboard::Down) {
                    visualizingScalingFactor = (visualizingScalingFactor <= 5) ? 5 : visualizingScalingFactor - 2;
                }

                if (visualizationWindowEvent.key.code == sf::Keyboard::V && state_solveSystem) {
                    std::cout << "Solving: " << attractors[attractorId] << std::endl;

                    // Identify the attractor to solve and visualize.
                    switch (attractorId)
                    {
                        case 0:
                            y0[0] = 1.1, y0[1] = 2.0, y0[2] = 7.0;    // initial value at time 0
                            tFinal = 50.0;
                            targetFunction = Attractor::Lorenz;
                            break;
                        case 1:
                            y0[0] = 1.1, y0[1] = 1.1, y0[2] = -0.01; // initial value at time 0
                            tFinal = 150.0;
                            targetFunction = Attractor::Thomas;
                            break;
                        case 2:
                            y0[0] = 0.1, y0[1] = 0.0, y0[2] = 0.0;  // initial value at time 0
                            tFinal = 50.0;
                            targetFunction = Attractor::Aizawa;
                            break;
                        case 3:
                            y0[0] = 1.1, y0[1]= 2.1, y0[2] = -2.0;  // initial value at time 0
                            tFinal = 50.0;
                            targetFunction = Attractor::Dadras;
                            break;
                        case 4:
                            y0[0] = 0.1, y0[1] = 0.0, y0[2] = 0.0; // initial value at time 0
                            tFinal = 100.0;
                            targetFunction = Attractor::Rossler;
                            break;
                        case 5:
                            y0[0] = 0.1, y0[1] = 0.0, y0[2] = 0.0; // initial value at time 0
                            tFinal = 120.0;
                            targetFunction = Attractor::SprottB;
                            break;
                        case 6:
                            y0[0] = 0.1, y0[1] = 0.0, y0[2] = 0.0; // initial value at time 0
                            tFinal = 50.0;
                            targetFunction = Attractor::Arneodo;
                            break;
                        case 7:
                            y0[0] = -0.2, y0[1] = -2.82, y0[2] = 2.71; // initial value at time 0
                            tFinal = 50.0;
                            targetFunction = Attractor::Lorenz83;
                            break;
                        case 8:
                            y0[0] = 0.1, y0[1] = 0.0, y0[2] = 0.0; // initial value at time 0
                            tFinal = 50.0;
                            targetFunction = Attractor::Halvorsen;
                            break;
                        case 9:
                            y0[0] = 0.1, y0[1] = 0.0, y0[2] = 0.0; // initial value at time 0
                            tFinal = 120.0;
                            targetFunction = Attractor::Sprott_LinzF;
                            break;
                    }

                    // Allocate the solver object
                    try // Sometimes operator new can fail.
                    {
                        std::cout << "Allocation successfull" << std::endl;
                        solver = new Solver(controller, targetFunction, y0, (size_t) ODE_DIMENSION, t0, tFinal, absTol, relTol, true);
                    }
                    catch(const std::bad_alloc &e)  // if new fails to allocate memory. Handle the error.
                    {
                        std::cerr << "Fail to allocate memory for the solver. " << e.what() << " Program close." << '\n';
                        visualizationWindow.close();
                    }

                    solver->solve();
                    std::cout << solver->m_yOut.size() << std::endl;

                    std::fill(limitPoints.begin(), limitPoints.end(), 0.0);

                    for (size_t i = 0; i < solver->m_yOut.size(); i++)
                    {
                        for (size_t j = 0; j < ODE_DIMENSION; j++)
                        {
                            long double point = std::fabsl(solver->m_yOut[i][j]);
                            if (point > limitPoints[j]) {
                                limitPoints[j] = point;
                            }
                        }
                    }

                    long double limitPoint = (limitPoints[targetAxis] < 0) ? std::floorl(limitPoints[targetAxis]) : std::ceill(limitPoints[targetAxis]);

                    axis1LabelText.setString(axis[targetAxis]);
                    tFinalLabelText.setString(std::string("Tf=" + std::to_string((int) tFinal)));
                    minPointLabelText.setString(std::string("-" + std::to_string(limitPoint)));
                    maxPointLabelText.setString(std::string("+" + std::to_string(limitPoint)));

                    float axis1LabelTextWidth = axis1LabelText.getLocalBounds().width;
                    float maxPointLabelTextHeight = maxPointLabelText.getLocalBounds().height;
                    axis1LabelText.setPosition(0.f, (PLOTTER_WINDOW_SIZE/2.f - (lineHeight - axis1LabelTextWidth / 2.f)));
                    maxPointLabelText.setPosition(sf::Vector2f(0.f, (PLOTTER_WINDOW_SIZE/2.f) - (2.f * lineHeight + 2.f * maxPointLabelTextHeight)));
                    
                    plottingScalingFactor = (plottingWindowSize / (float)tFinal);

                    state_solveSystem = false;
                    state_visualize = true;
                }
            }

            if (visualizationWindowEvent.type == sf::Event::MouseButtonPressed) {
                if (visualizationWindowEvent.mouseButton.button == sf::Mouse::Left) {
                    float xMPos = sf::Mouse::getPosition(visualizationWindow).x;
                    float yMPos = sf::Mouse::getPosition(visualizationWindow).y;

                    if ((xMPos >= attractorListBtnPos.x && xMPos <= attractorListBtnPos.x + attractorListBtnDim.x) && 
                        (yMPos >= attractorListBtnPos.y && yMPos <= attractorListBtnPos.y + attractorListBtnDim.y))
                    {
                        std::cout << "Button Clicked!" << std::endl;
                        state_isAttractorListBtnClicked = true;
                    }

                    if (state_isAttractorListBtnClicked) {
                        for (size_t i = 0; i < attractors.size(); i++)
                        {
                            float attractorListTextPosx = attractorListBtnPos.x;
                            float attractorListTextPosy = (attractorListBtnPos.y + attractorListBtnDim.y) + 16.f * i;

                            if ((xMPos >= attractorListTextPosx && xMPos <= attractorListTextPosx + attractorListBtnDim.x) &&
                                (yMPos >= attractorListTextPosy && yMPos <= attractorListTextPosy + 12.f))
                            {
                                std::cout << "You clicked: " << attractors[i] << std::endl;

                                // Deallocate the solver object
                                delete solver;       // deleting null pointer has no effect. Therefore checking is not needed here.
                                solver = nullptr;    // set the object to null pointer to avoid dangling pointer.

                                attractorId = i;
                                state_visualize = false;
                                state_solveSystem = true;
                                state_isAttractorListBtnClicked = false;
                                break;
                            }
                        }
                    }
                }
            }
        }

        sf::Event plottingWindowEvent;
        while (plottingWindow.pollEvent(plottingWindowEvent))
        {
            if (plottingWindowEvent.type == sf::Event::Closed)
            {
                plottingWindow.close();
            }

             if (plottingWindowEvent.type == sf::Event::MouseButtonPressed) {
                if (plottingWindowEvent.mouseButton.button == sf::Mouse::Left) {
                    float xMPos = sf::Mouse::getPosition(plottingWindow).x;
                    float yMPos = sf::Mouse::getPosition(plottingWindow).y;

                    if ((xMPos >= axisListBtnPos.x && xMPos <= axisListBtnPos.x + axisListBtnDim.x) && 
                        (yMPos >= axisListBtnPos.y && yMPos <= axisListBtnPos.y + axisListBtnDim.y))
                    {
                        std::cout << "Button Clicked!" << std::endl;
                        state_isAxisListBtnClicked = true;
                    }

                    if (state_isAxisListBtnClicked) {
                        for (size_t i = 0; i < ODE_DIMENSION; i++)
                        {
                            float axisListTextPosx = axisListBtnPos.x;
                            float axisListTextPosy = (axisListBtnPos.y + axisListBtnDim.y) + 16.f * i;

                            if ((xMPos >= axisListTextPosx && xMPos <= axisListTextPosx + axisListBtnDim.x) &&
                                (yMPos >= axisListTextPosy && yMPos <= axisListTextPosy + 12.f))
                            {
                                targetAxis = i;
                                state_isAxisListBtnClicked = false;
                                break;
                            }
                        }

                        long double limitPoint = (limitPoints[targetAxis] < 0) ? std::floorl(limitPoints[targetAxis]) : std::ceill(limitPoints[targetAxis]);
                        
                        axis1LabelText.setString(axis[targetAxis]);
                        minPointLabelText.setString(std::string("-" + std::to_string(limitPoint)));
                        maxPointLabelText.setString(std::string("+" + std::to_string(limitPoint)));
                        
                        float axis1LabelTextWidth = axis1LabelText.getLocalBounds().width;
                        float maxPointLabelTextHeight = maxPointLabelText.getLocalBounds().height;      
                        axis1LabelText.setPosition(0.f, (PLOTTER_WINDOW_SIZE/2.f - (lineHeight - axis1LabelTextWidth / 2.f)));
                        maxPointLabelText.setPosition(sf::Vector2f(0.f, (PLOTTER_WINDOW_SIZE/2.f) - (2.f * lineHeight + 2.f * maxPointLabelTextHeight)));
                    }
                }
            }

        }

        if (state_visualize) {
            std::vector<std::array<long double, 2>> projectedPoints(solver->m_yOut.size());
            // updateRotationMatrix(theta, xRotationMatrix, yRotationMatrix, zRotationMatrix);
            // doRotationAndProjection(solver->m_yOut, projectedPoints, orthogonalProjectionMatrix, xRotationMatrix, yRotationMatrix, zRotationMatrix);

            // Quaternion Rotation
            for (size_t i = 0; i < solver->m_yOut.size(); i++)
            {
                long double x = solver->m_yOut[i][0], y = solver->m_yOut[i][1], z = solver->m_yOut[i][2];
                Quaternion p(0.0, x, y, z);                                  // convert the point into pure quaternion
                Quaternion q = axisAngleToQuaternion(rotationAxis, theta);   // create a rotation quaternion
                Quaternion p_prime(0.0, 0.0, 0.0, 0.0);                      // stored the rotated points
                qRotation(p, q, p_prime);

                std::array<long double, 2> temp = {p_prime.x, p_prime.y};
                projectedPoints[i] = temp;
            }

            theta += 0.1;    // update theta

            // Draw on the screen
            visualizeODE(projectedPoints, visualizationWindow, visualizingScalingFactor);
            plotFunction(solver->m_yOut, solver->m_tOut, limitPoints, plottingWindow, plottingWindowOffsetX, plottingScalingFactor, lineHeight, targetAxis);

            plottingWindow.draw(axis1);
            plottingWindow.draw(axis2);
            plottingWindow.draw(tFinalLine);
            plottingWindow.draw(maxPointLine);
            plottingWindow.draw(axis1LabelText);
            plottingWindow.draw(axis2LabelText);
            plottingWindow.draw(originPointLine);
            plottingWindow.draw(tFinalLabelText);
            plottingWindow.draw(maxPointLabelText);
            plottingWindow.draw(minPointLabelText);
            plottingWindow.draw(originPointLabelText);
        }

        if (state_isAttractorListBtnClicked)
        {
            for (size_t i = 0; i < attractors.size(); i++)
            {
                if ((int) i == attractorId) {
                    attractorListText.setFillColor(sf::Color::Red);
                } else {
                    attractorListText.setFillColor(sf::Color::White);
                }

                attractorListText.setString(attractors[i]);
                attractorListTextPos.x = attractorListBtnPos.x;
                attractorListTextPos.y = (attractorListBtnPos.y + attractorListBtnDim.y) + 16 * i;
                attractorListText.setPosition(attractorListTextPos);
                visualizationWindow.draw(attractorListText);
            }
        }

        if (state_isAxisListBtnClicked)
        {
            for (size_t i = 0; i < ODE_DIMENSION; i++)
            {
                if ((int) i == targetAxis) {
                    axisListText.setFillColor(sf::Color::Black);
                } else {
                    axisListText.setFillColor(plottingWindowTextColor);
                }
                axisListText.setString(axis[i]);
                axisListTextPos.x = axisListBtnPos.x;
                axisListTextPos.y = (axisListBtnPos.y + axisListBtnDim.y) + 16 * i;
                axisListText.setPosition(axisListTextPos);
                plottingWindow.draw(axisListText);
            }
        }

        if (state_solveSystem) {
            std::string t = "Press 'v' to visualize " + attractors[attractorId] + " attractor."; 
            messageText1.setString(t);
            float messageText1Width = messageText1.getLocalBounds().width;
            float messageText1offSetX = (VISUALIZER_WINDOW_SIZE - messageText1Width) / 2.f;
            messageText1.setPosition(sf::Vector2f(messageText1offSetX, 100.f));
            visualizationWindow.draw(messageText1);
        }

        visualizationWindow.draw(attractorListBtnLabelText);
        visualizationWindow.draw(attractor_list_btn_sprite);
        visualizationWindow.draw(messageText2);

        plottingWindow.draw(axisListbtnLabelText);
        plottingWindow.draw(axis_list_btn_sprite);

        visualizationWindow.display();
        plottingWindow.display();
    }

    return 0;
}

void visualizeODE(std::vector<std::array<long double, 2>> &projectedPoints,
            sf::RenderWindow &targetWindow, float visualizingScalingFactor)
{
    for (size_t i = 0; i < projectedPoints.size()-1; i++) {
        sf::Vertex line[] =
        {
            sf::Vertex(sf::Vector2f(((VISUALIZER_WINDOW_SIZE/2) + projectedPoints[i][0] * visualizingScalingFactor), (VISUALIZER_WINDOW_SIZE/2) + (projectedPoints[i][1] * visualizingScalingFactor))),
            sf::Vertex(sf::Vector2f(((VISUALIZER_WINDOW_SIZE/2) + projectedPoints[i+1][0] * visualizingScalingFactor), (VISUALIZER_WINDOW_SIZE/2) + (projectedPoints[i+1][1] * visualizingScalingFactor)))
        };

        targetWindow.draw(line, 2, sf::Lines);
    }
}

void plotFunction(std::vector<std::vector<long double>> &y_out,
            std::vector<long double> &t_out, std::vector<long double> &limitPoints,
            sf::RenderWindow &targetWindow, float offsetX, float plottingScalingFactor, 
            float lineHeight, int targetAxis)
{

    for (size_t i = 0; i < y_out.size() - 1; i++)
    {
        /* code */
        long double limitPoint = limitPoints[targetAxis];
        limitPoint = (limitPoint < 0) ? std::floorl(limitPoint) : std::ceill(limitPoint);
        long double normVal1 = y_out[i][targetAxis] / limitPoint;
        long double normVal2 = y_out[i+1][targetAxis] / limitPoint;

        sf::Vertex line[] =
            {
                sf::Vertex(sf::Vector2f(offsetX + t_out[i] * plottingScalingFactor, (PLOTTER_WINDOW_SIZE/2 - lineHeight) + (normVal1 * lineHeight * -1.0))),
                sf::Vertex(sf::Vector2f(offsetX + t_out[i + 1] * plottingScalingFactor, (PLOTTER_WINDOW_SIZE/2 - lineHeight) + (normVal2 * lineHeight * -1.0)))
            };

        line[0].color = plottingLineColor;
        line[1].color = plottingLineColor;
        targetWindow.draw(line, 2, sf::Lines);
    }
}

long double vectorMag3(std::array<long double, ODE_DIMENSION> &vec)
{
    return std::sqrtl(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void normalize3D(std::array<long double, ODE_DIMENSION> &vec)
{
    long double mag = vectorMag3(vec);
    if (mag > 0.0) {
        vec[0] /= mag;
        vec[1] /= mag;
        vec[2] /= mag;
    }
}

void qRotation(Quaternion &p, Quaternion &q, Quaternion &p_prime)
{
    p_prime = q * p;
    q.conjugate();
    p_prime = p_prime * q;
}

Quaternion axisAngleToQuaternion(std::array<long double, ODE_DIMENSION> &axis, long double angleRadians)
{
    long double c = std::cosl(angleRadians/2.0);
    long double s = std::sinl(angleRadians/2.0);

    normalize3D(axis);

    return Quaternion(c, axis[0] * s, axis[1] * s, axis[2] * s);
}