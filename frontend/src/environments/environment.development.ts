export const environment = {
    production: false,
    firebase: {
        projectId: 'twocube-web',
        apiKey: 'REMOVED', // Not needed for emulator but required by Firebase SDK
        authDomain: 'app.twocube.ai',
        databaseURL: 'http://localhost:9000?ns=twocube-web', // Realtime DB (if applicable)
        storageBucket: 'twocube-web.appspot.com',
        messagingSenderId: 'your-messaging-sender-id',
        appId: 'your-app-id',
    
        // Use emulator ports
        emulator: {
        auth: 'http://localhost:9099',
        firestore: 'http://localhost:8080',
        functions: 'http://localhost:5001',
        hosting: 'http://localhost:5000',
        storage: 'http://localhost:9199',
        }
    },

    // API Endpoints - Replace with local development URLs
    chatAPIEndpoint: 'http://localhost:5001/twocube-web/us-central1/master_agent_interaction',
    l3chatAPIEndpoint: 'http://localhost:5001/twocube-web/us-central1/l3_master_agent_interaction',
    nameMakerAPIEndpoint: 'http://localhost:5001/twocube-web/us-central1/name_maker',

    // Box-related APIs
    boxRequestApi: 'http://localhost:5001/twocube-web/us-central1/request_sandbox',
    boxExecuteApi: 'http://localhost:5001/twocube-web/us-central1/execute_on_sandbox',
    boxCloseApi: 'http://localhost:5001/twocube-web/us-central1/close_sandbox',
    boxStatusApi: 'http://localhost:5001/twocube-web/us-central1/sandbox_status',
    boxUploadApi: 'http://localhost:5001/twocube-web/us-central1/upload_to_sandbox',
    addFirebaseFileApi: 'http://localhost:5001/twocube-web/us-central1/firebase_storage_to_sandbox',
    boxTerminalApi: 'http://localhost:5001/twocube-web/us-central1/run_terminal_command',
};
