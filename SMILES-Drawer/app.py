import http.server
import socketserver
import threading
import time
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

PORT = 8000

class SmilesRequestHandler(http.server.SimpleHTTPRequestHandler):
    """Handles GET requests for files and POST requests for SMILES data."""
    def do_POST(self):
        if self.path == '/save_smiles':
            try:
                content_length = int(self.headers['Content-Length'])
                post_data = self.rfile.read(content_length)
                smiles_string = post_data.decode('utf-8')
                
                print(f"[SMILES Received]: {smiles_string}")
                
                self.send_response(200)
                self.send_header('Content-type', 'text/plain')
                self.end_headers()
                self.wfile.write(b"OK")
            except Exception as e:
                print(f"Error processing POST request: {e}")
                self.send_response(500)
                self.end_headers()
                self.wfile.write(b"ERROR")
        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(b"Not Found")

    def do_GET(self):
        # Serve files using the default handler
        super().do_GET()

# --- Server Setup ---
HttdServer = socketserver.TCPServer(('', PORT), SmilesRequestHandler)

def run_server():
    """Runs the HTTP server in a separate thread."""
    print(f"Starting server on http://localhost:{PORT}")
    HttdServer.serve_forever()

# --- Main Application Logic ---
if __name__ == "__main__":
    # Start the server in a background thread
    # The daemon=True flag means the thread will automatically die when the main script exits
    server_thread = threading.Thread(target=run_server, daemon=True)
    server_thread.start()
    time.sleep(1) # Give the server a moment to start

    # --- Browser Setup ---
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    # This is the key option to create a minimal app-like window
    options.add_argument(f"--app=http://localhost:{PORT}")
    # Optional: define a window size
    options.add_argument("--window-size=1024,768")
    # Suppress Chrome error messages and warnings
    options.add_argument("--disable-logging")
    options.add_argument("--log-level=3")
    options.add_argument("--disable-web-security")
    options.add_argument("--disable-features=VizDisplayCompositor")
    options.add_argument("--disable-extensions")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_experimental_option('excludeSwitches', ['enable-logging'])
    options.add_experimental_option('useAutomationExtension', False)

    driver = None
    try:
        print("Launching application window...")
        print("Close the application window to shut down the server.")
        driver = webdriver.Chrome(service=service, options=options)
        
        # The script will now wait here until the browser window is closed.
        # We can do this by periodically checking if the window is still open.
        while True:
            if not driver.window_handles:
                # This means the user has closed the window
                break
            time.sleep(0.5)

    except Exception as e:
        # This block will catch errors if the user closes the window before the session starts
        # or if there's another WebDriver error.
        if "invalid session id" not in str(e) and "session deleted" not in str(e):
             print(f"An error occurred: {e}")

    finally:
        print("\nApplication window closed. Shutting down server...")
        if driver:
            driver.quit()
        HttdServer.shutdown() # Stop the server
        print("Server has been shut down. Goodbye.")