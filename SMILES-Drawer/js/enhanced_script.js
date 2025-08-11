var jsmeApplet; //Use var to make it a global variable accessible to selenium

function jsmeOnLoad() {
    //this function will be called after the JSME editor has been loaded
    jsmeApplet = new JSApplet.JSME("jsme_container", "800px", "500px", {
        "options" : "oldlook,star,reaction"
    });

    // Show status indicator
    const statusIndicator = document.getElementById("statusIndicator");
    const statusText = document.getElementById("statusText");
    
    function showStatus(message, type = 'info') {
        statusText.textContent = message;
        statusIndicator.className = `status-indicator status-${type}`;
        statusIndicator.style.display = 'block';
        
        // Auto-hide after 3 seconds unless it's a success message
        if (type !== 'success') {
            setTimeout(() => {
                statusIndicator.style.display = 'none';
            }, 3000);
        }
    }
    
    showStatus("Structure editor loaded successfully!", "success");

    // Get SMILES button
    document.getElementById("getSmilesBtn").onclick = function() {
        try {
            let smiles = jsmeApplet.smiles();
            document.getElementById("smilesOutput").value = smiles;
            
            if (smiles.trim()) {
                showStatus("SMILES generated successfully!", "success");
            } else {
                showStatus("No structure found. Please draw a molecule first.", "info");
            }
        } catch (error) {
            showStatus("Error generating SMILES", "error");
            console.error("Error getting SMILES:", error);
        }
    }

    // Save to Python button
    document.getElementById("saveSmilesBtn").onclick = function() {
        try {
            let smiles = jsmeApplet.smiles();
            
            if (!smiles.trim()) {
                showStatus("Please draw a structure first!", "info");
                return;
            }
            
            document.getElementById("smilesOutput").value = smiles; // Also update the text area
            showStatus("Sending structure to Python...", "info");

            // Send the SMILES string to the Python server
            fetch("/save_smiles", {
                method: "POST",
                headers: {
                    "Content-Type": "text/plain",
                },
                body: smiles,
            })
            .then(response => {
                if (response.ok) {
                    return response.text();
                } else {
                    throw new Error(`Server responded with status: ${response.status}`);
                }
            })
            .then(data => {
                console.log("Response from server:", data);
                showStatus("✓ Structure sent to Python successfully! You can now close this window.", "success");
                
                // Optional: Auto-close after a delay
                // setTimeout(() => {
                //     window.close();
                // }, 2000);
            })
            .catch((error) => {
                console.error("Error sending SMILES:", error);
                showStatus("Error: Could not connect to Python. Make sure the application is running.", "error");
            });
        } catch (error) {
            showStatus("Error processing structure", "error");
            console.error("Error in save function:", error);
        }
    }

    // Clear structure button
    document.getElementById("clearBtn").onclick = function() {
        try {
            jsmeApplet.clear();
            document.getElementById("smilesOutput").value = "";
            showStatus("Structure cleared", "info");
        } catch (error) {
            showStatus("Error clearing structure", "error");
            console.error("Error clearing:", error);
        }
    }

    // Add keyboard shortcuts
    document.addEventListener('keydown', function(event) {
        if (event.ctrlKey || event.metaKey) {
            switch(event.key) {
                case 'Enter':
                    event.preventDefault();
                    document.getElementById("getSmilesBtn").click();
                    break;
                case 's':
                    event.preventDefault();
                    document.getElementById("saveSmilesBtn").click();
                    break;
                case 'Delete':
                case 'Backspace':
                    if (event.shiftKey) {
                        event.preventDefault();
                        document.getElementById("clearBtn").click();
                    }
                    break;
            }
        }
    });
    
    // Show keyboard shortcuts hint
    console.log("Keyboard shortcuts:");
    console.log("Ctrl+Enter: Get SMILES");
    console.log("Ctrl+S: Save to Python");
    console.log("Ctrl+Shift+Delete: Clear structure");
}

// Error handling for JSME loading
window.addEventListener('error', function(event) {
    console.error('Script error:', event.error);
    const statusIndicator = document.getElementById("statusIndicator");
    const statusText = document.getElementById("statusText");
    
    if (statusIndicator && statusText) {
        statusText.textContent = "Error loading structure editor";
        statusIndicator.className = "status-indicator status-error";
        statusIndicator.style.display = 'block';
    }
});

// Add some CSS for error status
const style = document.createElement('style');
style.textContent = `
    .status-error {
        background-color: #ffe6e6;
        color: #d32f2f;
        border: 1px solid #ffcccc;
    }
`;
document.head.appendChild(style);
