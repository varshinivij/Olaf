# OLAF Backend & Firebase Cloud Functions

## Table of Contents

1. [Introduction](#introduction)  
2. [Features](#features)  
3. [Architecture Overview](#architecture-overview)  
   - [Agents, Routers, and Pipes](#agents-routers-and-pipes)  
   - [Data Flow Through the System](#data-flow-through-the-system)  
   - [Project Structure](#project-structure)  
4. [Getting Started](#getting-started)  
   - [Prerequisites](#prerequisites)  
   - [Installation](#installation)  
   - [Running with Docker (Recommended)](#running-with-docker-recommended)  
   - [Deployment](#deployment)  
   - [Usage](#usage)  
5. [MacOS Security Issues](#macos-security-issues)

---

## Introduction

This repository contains the code for Firebase Cloud Functions implemented in Python for the olaf project. Firebase Cloud Functions are used to execute backend code in response to events triggered by Firebase features and HTTPS requests.

This README demonstrates how to set up and deploy Python-based Firebase Cloud Functions for the olaf project. It also includes the engineering ideas behind our architecture as well as an illustrative flowchart, providing clarity on how agents, routers, and pipes integrate and how an API request flows through the system. The project folder structure is also outlined here.

---

## Features

- **Written in Python**  
- **Responds to HTTPS requests**  
- **Interacts with Firebase services**  
- **Scalable and serverless**  

---

## Architecture Overview

The backend architecture integrates multiple components to handle user requests, transform session data, and interact with AI agents. The primary elements are:

1. **Agents**: Produce responses from given input (prompt, history), possibly calling out to language models or tools.  
2. **Routers**: Orchestrate the flow of data by choosing which agent route to call, applying transformations via pipes, and chaining multiple routes if needed.  
3. **Pipes**: Modular transformations that preprocess or postprocess session data before it reaches an agent or after receiving an agent’s response.

### Agents, Routers, and Pipes

1. **Agent**  
   - Receives a prompt and conversation history.  
   - Interacts with tools (e.g., language models) and returns `(destination, response_generator)`.  
   - Focuses solely on producing a meaningful response.

2. **Pipe**  
   - Asynchronously applies transformations to session data (e.g., cleaning input, adding metadata, filtering sensitive content).  
   - Composable, so multiple pipes can be chained.  
   - Unaware of routing or agents; just transforms data.

3. **Router**  
   - Receives session data and a target route.  
   - Applies global pipes and route-specific pipes.  
   - Calls the assigned agent route function, retrieves `(destination, response_generator)`.  
   - If `destination == "user"`, returns the response; otherwise, re-runs the process on the specified route.

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

### Project Structure

Below is a brief overview describing the structure and purpose of each folder, loosely based on Angular:

- `/agents`: Classes representing agents, to be used by routes.  
- `/datastructures`: Classes containing data structures, including the Router implementation.  
- `/endpoints`: Modules containing all public-facing endpoints exposed to Firebase Functions.  
- `/models`: Classes modeling database entities.  
- `/pipers`: Classes representing pipes, to be used in the Router. (`/pipes` is reserved by Firebase.)  
- `/routes`: Modules containing routes as functions, each using agents, returning responses as generators, and specifying a destination.  
- `/services`: Classes containing utility methods, database calls, and API calls for endpoints.  
- `/utils`: Miscellaneous utility files.  

Most endpoints use their respective services, while `agent_endpoints.py` streams responses using the Router data structure and routes imported from `/routes`.

---

## Getting Started

Follow these instructions to set up and deploy the Firebase Cloud Functions for your own OLAF project

### Prerequisites

- Node.js (to use the Firebase CLI)  
- Firebase CLI  
- Python 3.7 or later  
- A Firebase project  

### Installation

1. **Clone the repository:**
   ```sh
   git clone [url to Olaf repo]
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
   # during initialization, select Python as your language and follow the prompts
   ```

5. **Set up a virtual environment and install dependencies:**
   ```sh
   python -m venv venv
   source venv/bin/activate  # on Windows use `venv\Scripts\activate`
   pip install -r requirements.txt
   ```

### Running with Docker (Recommended)

To simplify local development, you can use Docker to run the Firebase Emulator with your cloud functions.

1. **Build the Docker image:**
   ```sh
   docker compose build
   ```

2. **Run the container:**
   ```sh
   docker compose up
   ```

3. **Check logs:**
   ```sh
   docker compose logs -f
   ```

4. **Stop the emulator:**
   ```sh
   docker compose down
   ```

5. **Alternative Easy Scripts**
   ```sh
   /.run_firebase
   ```
   ```sh
   ./stop_firebase
   ```

### Deployment

1. **Set the working directory:**
   ```sh
   cd functions
   ```

2. **Deploy the functions to Firebase:**
   ```sh
   firebase deploy --only functions
   ```

3. **OR deploy to the Firebase emulator for testing:**
   ```sh
   firebase emulators:start --only functions
   ```

### Usage

To test and use the deployed functions, follow these steps:

1. **Trigger a function via HTTPS:**  
   Use tools like `curl`, Postman, or a web browser to send a request to your function’s URL. The URL can be found in the Firebase console under the Functions section.  
   ```sh
   curl -X GET https://us-central1-your-project-id.cloudfunctions.net/yourFunctionName
   ```

2. **Event-driven functions:**  
   Set up Firebase events (such as Firestore triggers) to invoke the functions automatically based on database changes.

---

## MacOS Security Issues

On macOS, you may need to run:

```sh
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
```

This environment variable allows you to serve functions locally over the emulator on macOS environments.