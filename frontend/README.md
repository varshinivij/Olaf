# OLAF (Open Life Science Analysis Framework) Frontend

Welcome to **OLAF (Open Life Science Analysis Framework) Frontend**, an open-source Angular project! 🚀 This repository is designed to streamline bioinformatics workflows with an intuitive frontend interface. Contributions are always welcome! ❤️

## Table of Contents

- [About the Project](#about-the-project)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Development Server](#development-server)
- [Project Structure](#project-structure)
- [Building the Project](#building-the-project)
- [Testing](#testing)
  - [Running Unit Tests](#running-unit-tests)
  - [Running End-to-End Tests](#running-end-to-end-tests)
---

## About the Project

The OLAF frontend is built with Angular and provides a robust UI for LLM driven Bioinformatics applications. It is designed to be modular, scalable, and easily extendable.

## Getting Started

### Prerequisites

Ensure you have the following installed:

- [Node.js](https://nodejs.org/) (latest LTS version recommended)
- [Angular CLI](https://angular.dev/tools/cli)

### Installation

Clone the repository and install dependencies:

```sh
# Clone the repository
git [olaf repo url]
cd olaf/frontend

# Install dependencies
npm install
```

### Development Server

Run the following command to start a local development server:

```sh
ng serve
```

Navigate to `http://localhost:4200/` in your browser. The application will automatically reload when you modify source files.

## Project Structure

```
Olaf/
├── src/
│   ├── app/                # Main application code
│   ├── environments/       # Environment configurations
│   ├── assets/             # Static assets
│   ├── styles/             # Global styles
│   ├── main.ts             # Application entry point
│   └── index.html          # Main HTML file
└── angular.json            # Angular configuration file
```

## Building the Project

To build the project, run:

```sh
ng build
```

The build artifacts will be stored in the `dist/` directory.

## Testing

### Running Unit Tests

Execute the unit tests using [Karma](https://karma-runner.github.io/):

```sh
ng test
```

### Running End-to-End Tests

To run end-to-end tests:

```sh
ng e2e
```

Ensure you have a platform-specific e2e testing package installed before running this command.

