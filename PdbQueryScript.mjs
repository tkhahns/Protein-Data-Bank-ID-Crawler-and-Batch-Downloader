import fetch from 'node-fetch';

const baseUrl = 'https://search.rcsb.org/rcsbsearch/v2/query';

const data = encodeURIComponent(
    JSON.stringify({
        query: {
            type: "terminal",
            service: "text"
            },
            return_type: "entry"
    }))
const url = baseUrl + `?json=${data}`

fetch(url)
    .then(res => res.text())
    .then(text => console.log(text));
