# TwoCubeFrontend

Welcome to **TwoCubeFrontend**, an open-source Angular project! 🚀 This repository is designed to streamline bioinformatics workflows with an intuitive frontend interface. Contributions are always welcome! ❤️

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
- [Firebase Setup](#firebase-setup)
- [Requesting Access to Google Cloud](#requesting-access-to-google-cloud)
- [Contributing](#contributing)
- [License](#license)

---

## About the Project

TwoCubeFrontend is built with Angular and provides a robust UI for LLM driven Bioinformatics applications. It is designed to be modular, scalable, and easily extendable.

## Getting Started

### Prerequisites
Ensure you have the following installed:
- [Node.js](https://nodejs.org/) (latest LTS version recommended)
- [Angular CLI](https://angular.dev/tools/cli)

### Installation
Clone the repository and install dependencies:

```sh
# Clone the repository
git clone https://github.com/TwoCubeAI/TwoCubeFrontend.git
cd TwoCubeFrontend

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
TwoCubeFrontend/
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

## Firebase Setup

### 1. Configure Firebase Credentials
For security, Firebase credentials should be stored in a local `environment.secret.ts` file, which is excluded from version control.

#### Step 1: Create the file
Create a new file at `src/environments/environment.secret.ts`

#### Step 2: Add Firebase credentials
```js
export const secrets = {
    production: false,
    firebase: {
        projectId: 'YOUR_PROJECT_ID',
        appId: 'YOUR_APP_ID',
        storageBucket: 'YOUR_STORAGE_BUCKET',
        apiKey: process.env.FIREBASE_API_KEY || 'YOUR_FIREBASE_API_KEY',
        authDomain: 'YOUR_AUTH_DOMAIN',
        messagingSenderId: 'YOUR_MESSAGING_SENDER_ID',
        measurementId: 'YOUR_MEASUREMENT_ID',
    }
};
```

#### Step 3: Add to `.gitignore`
Ensure that `.gitignore` includes:

```
src/environments/environment.secret.ts
```

## Requesting Access to Google Cloud
If you are part of the core team and need access to the production Firebase environment, contact **Dylan Riffle**.

Access includes:
- Production Firebase credentials
- Firestore database permissions
- Cloud Functions deployment rights
- Access to Firebase Hosting and related services

📌 Note: Only verified team members will be granted access.

## Contributing
We ❤️ open-source contributions! If you'd like to contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature-name`)
3. Commit your changes (`git commit -m "Add feature"`)
4. Push to your fork (`git push origin feature-name`)
5. Open a pull request 🚀

For major changes, please open an issue first to discuss your proposal.

## License
This project is licensed under the **MIT License**
---

🚀 **Happy coding, and welcome to the TwoCubeFrontend community!** 🎉

