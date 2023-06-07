# Nonlinear Control
## A Repository for exploring topics in nonlinear control on simple systems

This repository will examine various topics in nonlinear control, with a focus on concepts presented in CDS 233 (Nonlinear Control) as taught by Dr. Aaron Ames at Caltech. It will focus on implementing controllers on simple, easy to compute with systems, while highlighting strengths and weaknesses of various control techniques. The repository is divided into several folders, each of which investigate different aspects of nonlinear control.

## Models
The Models folder contains the system dynamics for the inverted pendulum and cartpole, the two systems which will be considered in this Repository. The pendulum will be uses as the canonical example of a fully actuated nonlinear system, while the cartpole will serve to demonstrate the effects of underactuation. 

## QPs
The QPs folder contains the generic implementation of several of the quadratic program based controllers presented CDS 233, which enforce Control Lyapunov Function or Control Safety Function constraints to guarantee stability or safety of the system. It also contains, when possible, explicit solutions to these QPs which can be implemented without using an optimization program.

## Simulators
The simulators folder contains wrapper functions which simulate systems from the models folder with controllers provided by the user (potentially a controller wrapping one of the QPs from the QPs folder). The purpose of these simulation wrappers is to provide easy access to various implementation styles and uncertainties, including:
* Model Uncertainty
* Sample and Hold Inputs
* Input Delays (TODO)
* Sensor Noise (TODO)
* Disturbances (TODO)

These simulators aim to help evaluate the robustness of the given controllers to the uncertainties present in real-world implementation. Controllers which leverage nonlinear dynamics are inherently model-based, and as such care must be taken to ensure robustness to uncertainty and noise.

## Nominal Control
The nominal control folder contains methods for controlling nominal systems. It is primarily concerned with situations where the model is well known, state and input bounds are not a major concern, and sensor noise and disturbances are negligible (i.e. this situation is never the case). It focuses on comparing the results of linear (PD) controllers, feedback linearization controllers, control Lyapunov function based controllers, and (TODO) model predictive controllers.
