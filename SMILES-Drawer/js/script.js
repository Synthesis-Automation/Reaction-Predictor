var jsmeApplet; //Use var to make it a global variable accessible to selenium

function jsmeOnLoad() {
    //this function will be called after the JSME editor has been loaded
    jsmeApplet = new JSApplet.JSME("jsme_container", "800px", "500px", {
        "options" : "oldlook,star,reaction"
    });

    document.getElementById("getSmilesBtn").onclick = function() {
        let smiles = jsmeApplet.smiles();
        document.getElementById("smilesOutput").value = smiles;
    }

    document.getElementById("saveSmilesBtn").onclick = function() {
        let smiles = jsmeApplet.smiles();
        document.getElementById("smilesOutput").value = smiles; // Also update the text area

        // Send the SMILES string to the Python server
        fetch("/save_smiles", {
            method: "POST",
            headers: {
                "Content-Type": "text/plain",
            },
            body: smiles,
        })
        .then(response => response.text())
        .then(data => {
            console.log("Response from server:", data);
            // You could add a visual confirmation here if needed
        })
        .catch((error) => {
            console.error("Error sending SMILES:", error);
            alert("Error: Could not connect to the Python server. Is it running?");
        });
    }
}
