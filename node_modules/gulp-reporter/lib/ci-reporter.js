'use strict';
const checkstyleFormatter = require('checkstyle-formatter');
const reportBuilder = require('junit-report-builder');
const yaml = require('js-yaml');
const fs = require('fs-extra');
const axios = require('axios');
const ci = require('ci-info');
const path = require('path');
const url = require('url');
const os = require('os');
const icons = require('./icons');

const privatePrefix = ci.isCI && ci.name && ci.name.replace(/\s*CI$/i, '').toUpperCase() + '_';
const reportPath = path.join.bind(
	path,
	getEnv('REPORTS', 'REPORT_PATH', 'TEST_REPORTS') || circleReportPath() || 'lint-reports'
);

let appveyorApiUrl;

const category = {
	warn: 'Warning',
	info: 'Information',
	error: 'Error',
};
const severity = {
	warn: 'warning',
	info: 'information',
};

function noop () {
	return Promise.resolve();
}

function appveyor (errors) {
	return Promise.all(
		errors.map(error => (
			axios.post(appveyorApiUrl, {
				message: `${icons[error.severity] || icons.success} ${error.message}`,
				category: category[error.severity] || error.severity,
				details: [error.plugin, error.rule, (error.docShort || error.doc)].filter(Boolean).join(' '),
				fileName: error.fileName,
				line: error.lineNumber,
				column: error.columnNumber,
				projectName: [error.plugin, error.rule].filter(Boolean).join(' '),
			})
		))
	);
}

function getEnv (names) {
	names = Array.from(arguments);
	let name;
	let value;
	while ((name = names.pop())) {
		value = (privatePrefix && process.env[privatePrefix + name]) || process.env['CI_' + name] || process.env[name];
		if (value) {
			return resolveHomePath(value);
		}
	}
}

function junit (errors) {
	const builder = reportBuilder.newBuilder();
	Object.keys(errors).forEach(fileName => {
		const suite = builder.testSuite().name(fileName);
		errors[fileName].forEach(error => {
			const testCase = suite.testCase()
				.name(`${icons[error.severity] || icons.success} ${error.message}`)
				.className([error.plugin, error.rule].filter(Boolean).join(' '));
			const message = `${error.fileName}:${error.lineNumber}:${error.columnNumber}`;
			if (!error.severity || error.severity === 'error') {
				testCase.error(message);
			} else {
				testCase.failure(message);
			}
		});
	});

	return fs.outputFile(reportPath(String(Date.now()), 'lint-result.xml'), builder.build());
}

function checkstyle (errors) {
	const results = [];

	Object.keys(errors).forEach(filename => {
		results.push({
			filename,
			messages: errors[filename].map(error => ({
				severity: severity[error.severity] || error.severity,
				column: error.columnNumber,
				line: error.lineNumber,
				message: `${icons[error.severity] || icons.success} ${error.message}`,
				source: [error.plugin, error.rule].filter(Boolean).join(' '),
			})),
		});
	});

	return fs.outputFile(
		reportPath(String(Date.now()), 'checkstyle-result.xml'),
		checkstyleFormatter(results)
	);
}

function getErrorsByFile (errors) {
	const result = {};
	errors.forEach((error) => {
		if (result[error.fileName]) {
			result[error.fileName].push(error);
		} else {
			result[error.fileName] = [error];
		}
	});
	return result;
}

function reduceErrors (files) {
	return files.reduce((errors, file) => (
		errors.concat(file.report.errors)
	), []).filter(Boolean);
}

if (ci.isCI && (ci.APPVEYOR || ci.CIRCLE || ci.JENKINS)) {
	appveyorApiUrl = process.env.APPVEYOR_API_URL;
	if (appveyorApiUrl) {
		// eslint-disable-next-line node/no-deprecated-api
		appveyorApiUrl = url.resolve(appveyorApiUrl, 'api/build/compilationmessages');
	}
	module.exports = function (files) {
		const errors = reduceErrors(files);
		if (!errors.length) {
			return noop();
		}
		let result;
		if (appveyorApiUrl) {
			result = appveyor(errors);
		} else {
			const errorSet = getErrorsByFile(errors);
			if (ci.CIRCLE) {
				result = junit(errorSet);
			} else {
				result = checkstyle(errorSet);
			}
		}
		return result.then(noop);
	};
} else {
	module.exports = noop;
}

function resolveHomePath (strPath) {
	return strPath.replace(
		/^~(?=\/|$)/,
		os.homedir
	);
}

// https://circleci.com/docs/2.0/collect-test-data/
function circleReportPath () {
	const CIRCLE_JOB = process.env.CIRCLE_JOB;
	const CIRCLE_WORKING_DIRECTORY = process.env.CIRCLE_WORKING_DIRECTORY;
	if (!(ci.CIRCLE && CIRCLE_JOB && CIRCLE_WORKING_DIRECTORY)) {
		return;
	}
	const cwd = resolveHomePath(CIRCLE_WORKING_DIRECTORY);
	let config = path.join(
		cwd,
		'.circleci/config.yml'
	);
	config = fs.readFileSync(
		config,
		'utf-8'
	);
	config = yaml.safeLoad(config);

	config = config.jobs[CIRCLE_JOB].steps.reduce((results, step) => {
		step = step.store_test_results;
		if (step) {
			step = step.paths || step.path;
			if (step) {
				results = results.concat(step);
			}
		}
		return results;
	}, []);

	config = config.find(path => !/\w+\.\w+$/.test(path)) || config[0];

	if (config) {
		return path.resolve(
			cwd,
			resolveHomePath(config)
		);
	}
}
