import fetch from 'node-fetch';
import fs from 'fs';

const numOfStructures = 215908;
const pageLimit = 10000;
const baseUrl = 'https://search.rcsb.org/rcsbsearch/v2/query';

let ids = [];
let data = {
    query: {
        type: "terminal",
        service: "text"
    },
    request_options: {
        paginate: {
            start: 0,
            rows: pageLimit
        }
    },
    return_type: "entry"
};

const iterations = Math.floor(numOfStructures/pageLimit);
for (let i = 0; i < iterations; i++) {
    data.request_options.paginate.start = i * pageLimit;
    const url = baseUrl + `?json=${encodeURIComponent(JSON.stringify(data))}`;

    let currentIds = await fetch(url)
        .then(res => res.json())
        .then(json => json["result_set"].map(protein => protein["identifier"]));
    
    ids = ids.concat(currentIds);
}

const outputPath = 'list_file.txt';
fs.writeFile(outputPath, String(ids), err => {
    if (err) {
        console.error(err);
    }
})
