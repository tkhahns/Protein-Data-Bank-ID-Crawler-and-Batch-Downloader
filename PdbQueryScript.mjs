import fetch from 'node-fetch';

const numOfStructures = 215908;
const pageLimit = 10000;
const numOfPages = Math.floor(215908/100);
const baseUrl = 'https://search.rcsb.org/rcsbsearch/v2/query';

//let startNumbers = Array.from({length: numOfPages}, (_, i) => i * pageLimit);
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

// correct the number of loops
const iterations = Math.floor(numOfStructures/pageLimit);
for (let i = 0; i < 2; i++) {
    data.request_options.paginate.start = i * pageLimit;
    const url = baseUrl + `?json=${encodeURIComponent(JSON.stringify(data))}`;

    let currentIds = await fetch(url)
        .then(res => res.json())
        .then(json => json["result_set"].map(protein => protein["identifier"]));
    
    ids = ids.concat(currentIds);
}
console.log(ids);