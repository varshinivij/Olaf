Below is an updated README that includes the engineering ideas discussed, as well as the illustrative flowchart. This added section provides clarity on how agents, routers, and pipes integrate, and how an API request flows through the system.

# Two Cube Firebase Cloud Function

This repository contains the code for Firebase Cloud Functions implemented in Python for the Two Cube project. Firebase Cloud Functions are used to execute backend code in response to events triggered by Firebase features and HTTPS requests.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Deployment](#deployment)
- [Usage](#usage)
- [Architecture Overview](#architecture-overview)
  - [Agents, Routers, and Pipes](#agents-routers-and-pipes)
  - [Data Flow Through the System](#data-flow-through-the-system)
- [Mac Security Issues](#mac-security-issues)

## Introduction

Firebase Cloud Functions allow you to run server-side code in response to events triggered by Firebase features and HTTPS requests. This repository demonstrates how to set up and deploy Python-based Firebase Cloud Functions for the Two Cube project.

## Features

- Written in Python
- Responds to HTTPS requests
- Interacts with Firebase services
- Scalable and serverless

## Getting Started

Follow these instructions to set up and deploy the Firebase Cloud Functions for your Two Cube project.

### Prerequisites

- Node.js (to use the Firebase CLI)
- Firebase CLI
- Python 3.7 or later
- A Firebase project

### Installation

1. **Clone the repository:**

   ```sh
   git clone https://github.com/yourusername/two-cube-firebase-cloud-function.git
   cd two-cube-firebase-cloud-function
   ```
2. **Install Firebase CLI:**

   ```sh
   npm install -g firebase-tools
   ```
3. **Log in to Firebase:**

   ```sh
   firebase login
   ```
4. **Initialize Firebase in your project directory:**

   ```sh
   firebase init functions
   ```

   During initialization, select Python as your language and follow the prompts.
5. **Set up a virtual environment and install dependencies:**

   ```sh
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   pip install -r requirements.txt
   ```

### Deployment

1. **Deploy the functions to Firebase:**

   ```sh
   firebase deploy --only functions
   ```

2. **Run the Firebase emulator:**

   ```sh
   firebase emulators:start --only functions
   ```

### Usage

To test and use the deployed functions, follow these steps:

1. **Trigger a function via HTTPS:**

   Use tools like `curl`, Postman, or a web browser to send a request to your function's URL. The URL can be found in the Firebase console under the Functions section.

   ```sh
   curl -X GET https://us-central1-your-project-id.cloudfunctions.net/yourFunctionName
   ```
2. **Event-driven functions:**

   Set up Firebase events (such as Firestore triggers) to invoke the functions automatically based on database changes.

## Architecture Overview

In this backend architecture, we integrate multiple components to handle user requests, transform session data, and interact with AI agents. The primary elements are:

- **Agents**: Produce responses from a given input (prompt, history), possibly calling out to language models or tools.
- **Routers**: Orchestrate the flow of data by choosing which agent route to call, applying transformations via pipes, and chaining multiple routes if needed.
- **Pipes**: Modular transformations that preprocess or postprocess session data before it reaches an agent or after receiving an agent’s response.

### Agents, Routers, and Pipes

1. **Agent**:
   - Receives a prompt and conversation history.
   - Interacts with tools (e.g., language models) and returns `(destination, response_generator)` results.
   - Doesn't worry about routing or where data came from, just focuses on producing a meaningful response.
2. **Pipe**:
   - Applies transformations to session data (e.g., cleaning input, adding metadata, filtering sensitive content).
   - Is asynchronous and composable, so multiple pipes can be chained.
   - Doesn’t know about routes or agents; just transforms data.
3. **Router**:
   - Receives session data and a target route.
   - Applies global pipes and route-specific pipes.
   - Calls the assigned agent route function, retrieves `(destination, response_generator)`.
   - If `destination == "user"`, returns the response to the caller.
   - If `destination` is another route, re-runs the process to chain multiple agents/pipes.

By decoupling these components, new features or data flows can be introduced without rewriting the entire system. Developers can add or remove pipes and agents easily, and routes can be reconfigured as needed.

### Data Flow Through the System

Below is a conceptual flowchart showing how a request flows from an API request to the user response, incorporating agents, routers, and pipes:

```python
[HTTP Request] ---> [API Endpoint] ---> [session_data] ---> [Router]
    |                                                           |
    |                                                           v
    |                                              [Pipes (global, route-specific)]
    |                                                           |
    |                                                           v
    |                                                   [MasterAgent Route]
    |                                                           |
    |                                                           v
    |                                               [MasterAgent] ---> (destination, generator)
    |                                                   if destination = 'user':
    |                                                       return generator to API
    |                                                   if destination = another_route:
    |                                                       router.route(another_route, session_data)
    |                                                       *process continues*
    |
    v
[HTTP Response to User]
```

This architecture ensures that:

- Agents focus on generating responses.
- Pipes focus on data transformations.
- The Router focuses on orchestration, making it straightforward to trace the flow of data for new developers.

## Mac Security Issues

On macOS, you may need to run:

```sh
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
```

This environment variable allows you to serve functions locally over the emulator on macOS environments.
