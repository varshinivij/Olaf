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
- [Contributing](#contributing)
- [License](#license)

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

### Usage

To test and use the deployed functions, follow these steps:

1. **Trigger a function via HTTPS:**

    Use tools like `curl`, Postman, or a web browser to send a request to your function's URL. The URL can be found in the Firebase console under the Functions section.

    ```sh
    curl -X GET https://us-central1-your-project-id.cloudfunctions.net/yourFunctionName
    ```

2. **Event-driven functions:**

    Set up Firebase events (such as Firestore triggers) to invoke the functions automatically based on database changes.

### Mac Security issues

Some users may need to run:

    ```sh
    export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
    ```

For their macs to allow them to serve functions locally over the emulator.

### Example Function

Here is a basic example of a Firebase Cloud Function written in Python:

```python
from firebase_functions import https_fn

@https_fn.on_request
def hello_world(request):
    return https_fn.Response('Hello, World!', status=200)