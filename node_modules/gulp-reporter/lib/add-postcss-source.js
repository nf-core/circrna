'use strict';
const get = require('lodash.get');

function addPostcssSource (file, errors) {
	const source = get(file, 'postcss.root.source.input.css');
	if (source && errors.length) {
		const doc = source && source.match(/^.*$/gm);
		errors = errors.map(error => {
			const source = error.line && doc[error.line - 1];
			if (source) {
				error = Object.assign({
					source,
				}, error);
			}
			return error;
		});
	}
	return errors;
}

module.exports = addPostcssSource;
