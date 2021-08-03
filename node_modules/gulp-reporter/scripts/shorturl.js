'use strict';
// process.env.HTTP_PROXY = 'http://127.0.0.1:1080/';
const JSDOM = require('jsdom').JSDOM;
const stringify = require('json-stable-stringify');
const fs = require('fs-extra');
const got = require('got');
const googl = require('goo.gl');
const pkg = require('../package.json');
const isCI = require('ci-info').isCI;
const http = require('http');
const https = require('https');
// Set a developer key (_required by Google_; see http://goo.gl/4DvFk for more info.)
googl.setKey('AIzaSyACqNSi3cybDvDfWMaPyXZEzQ6IeaPehLE');

function awaitArray (arr) {
	return Promise.all(arr.map(item => {
		if (Array.isArray(item)) {
			return awaitArray(item);
		} else {
			return item;
		}
	}));
}

function expandArray (arr, result) {
	arr.forEach(item => {
		if (Array.isArray(item)) {
			expandArray(item, result);
		} else if (item) {
			result.push(item);
		}
	});
	return result;
}

function shortUrl (url) {
	return googl.shorten(url);
}

function shortUrlCn (url) {
	return got(`https://api.t.sina.com.cn/short_url/shorten.json?source=3271760578&url_long=${encodeURIComponent(url)}`, {
		json: true,
	}).then(result => {
		return result.body[0].url_short.replace(/^https?/i, 'https');
	});
}

const chechkLinkCache = {};
function chechkLinkWithCache (url) {
	if (!chechkLinkCache[url]) {
		chechkLinkCache[url] = chechkLink(url);
	}
	return chechkLinkCache[url];
}
function chechkLink (url) {
	if (/^https?:\/\/goo.gl\//i.test(url)) {
		// if (isCI) {
		return googl.expand(url).then(chechkLinkWithCache);
		// } else {
		// return Promise.reject(new Error(''));
		// }
	}
	if (url.startsWith('https://stylelint.io/') || /^https:\/\/(?:cn\.)?eslint\.org\//i.test(url)) {
		return Promise.resolve(url);
	}
	return new Promise((resolve, reject) => {
		function reTry () {
			req.abort();
			chechkLink(url).then(resolve, reject);
		}
		const req = (
			/^https:/i.test(url) ? https : http
		).get(url, res => {
			req.abort();
			if (res.statusCode >= 400) {
				reject(res.statusCode);
			} else if (res.statusCode >= 300) {
				if (res.headers.location) {
					chechkLinkWithCache(res.headers.location).then(resolve, reject);
				} else {
					reject(res.statusCode);
				}
			} else if (res.statusCode >= 200) {
				resolve(url);
			}
			res.resume();
		}).on(
			'error',
			reTry
		).setTimeout(
			0x1FFFF,
			reTry
		);
	});
}

let count = 0;

function get (url, selector) {
	return got(
		url
	).then(
		res => new JSDOM(res.body, {
			url: url,
			contentType: res.headers['content-type'],
			referrer: res.headers['referer'],
		}).window
	).then(window => {
		if (!isCI) {
			console.log('ready:', ++count, window.document.location.href);
		}
		return Array.from(window.document.querySelectorAll(selector)).map(a => a.href);
	}).catch(error => {
		if (!url.endsWith('/tree/HEAD/docs/rules')) {
			console.error(error);
		}
		return [];
	});
}

function updateFile (file, urls, shortUrlFn) {
	file = require.resolve(file);
	let hasChange = false;
	return fs.readJSON(file).then(shorturlCache => (
		Promise.all(
			urls.filter(
				(url, i) => !(url.toLowerCase() in shorturlCache)
			).filter(
				(url, i) => i < 100
			).map(
				url => chechkLinkWithCache(url).then(() => (
					shortUrlFn(url).then(shortUrl => {
						if (shortUrl) {
							hasChange = true;
							shorturlCache[url.toLowerCase()] = shortUrl;
						}
					}).catch(console.error)
				)).catch(() => {})
			)
		).then(() => (
			hasChange && fs.writeFile(
				file,
				stringify(
					shorturlCache,
					{
						space: '\t',
					}
				) + '\n',
				'utf8'
			)
		)).then(() => hasChange)
	));
}

const eslintRules = Object.keys(
	require('eslint/lib/load-rules')()
).map(
	rule => rule.toLowerCase()
);

const npmUrlPrefix = 'https://www.npmjs.com/package/eslint-plugin-';

const EslintPluginDocBaseUrl = {
	sql: rule => npmUrlPrefix + `sql#eslint-plugin-sql-rules-${rule}`,
	standard: () => npmUrlPrefix + 'standard#rules-explanations',
	gettext: rule => npmUrlPrefix + `gettext#gettext${rule}`,
};
[
	'flowtype',
	'jsdoc',
	'alint',
].forEach((plugin) => {
	EslintPluginDocBaseUrl[plugin] = rule => npmUrlPrefix + plugin + '#' + rule;
});

function getJSON (url) {
	return got(url, {
		json: true,
	}).then(
		res => res.body
	);
}

Promise.all([
	isCI && get('https://github.com/editorconfig/editorconfig/wiki/EditorConfig-Properties', '.markdown-body h3 a[href^="#"]'),
	// isCI && get('https://cn.eslint.org/docs/rules/', '.rule-list a[href]'),
	// isCI && get('https://eslint.org/docs/rules/', '.rule-list a[href]'),
	// isCI && get('https://stylelint.io/user-guide/rules/', 'h1 ~ ul a[href$="/"]'),
	isCI && get('https://palantir.github.io/tslint/rules/', '.rules-list a[href]').then(urls => (
		urls.map(url => (
			url.replace(/\/*$/, '/')
		))
	)),
	isCI && get('https://github.com/yaniswang/HTMLHint/wiki/Rules', '.markdown-body ul a[href*="HTMLHint"]'),
	isCI && get('https://github.com/CSSLint/csslint/wiki/Rules', '.markdown-body ul a[href^="/CSSLint"]'),
	isCI && get('http://jscs.info/rules', '.rule-list a[href]').then(urls => (
		urls.map(
			url => url.replace(/#/, '/')
		)
	)),

	// ESLint (zh-CN)
	eslintRules.map(rule => (
		`https://cn.eslint.org/docs/rules/${rule}`
	)),

	// ESLint
	eslintRules.map(rule => (
		`https://eslint.org/docs/rules/${rule}`
	)),

	// eslint-plugin-standard
	npmUrlPrefix + 'standard#rules-explanations',

	// eslint-plugin-sql
	npmUrlPrefix + 'sql#eslint-plugin-sql-rules-format',
	npmUrlPrefix + 'sql#eslint-plugin-sql-rules-no-unsafe-query',

	// eslint-plugin-compat
	[
		'serviceworker',
		'intersectionobserver',
		'webassembly',
		'paymentrequest',
		'serviceworker',
		'fetch',
		'promise',
	].map(s => 'https://www.caniuse.com/#search=' + s),

	// ESLint plugins docs from local node modules
	Object.keys(pkg.devDependencies).filter(
		pkgName => /^eslint-plugin-/.test(pkgName)
	).map(pkgName => {
		let baseUrl = EslintPluginDocBaseUrl[pkgName.slice(14)];
		if (!baseUrl) {
			let repository = require(pkgName + '/package.json').repository;
			repository = repository.url || repository;
			repository = repository.replace(/^(?:git\+)?\w+:\/+(?:.+?@)?(.+?)(?:\/+|\.git)?$/i, 'https://$1');
			baseUrl = rule => `${repository}/blob/HEAD/docs/rules/${rule}.md#readme`;
		}
		return Object.keys(require(pkgName).rules).map(baseUrl);
	}),

	// ESLint plugins docs from registry.npmjs.com
	isCI && getJSON('http://registry.npmjs.com/-/v1/search?text=eslint-plugin-&size=250').then(
		results => results.objects.map(
			object => object.package
		).filter(
			pkgData => pkgData.name.startsWith('eslint-plugin-')
		).filter(
			pkgData => !/^eslint-plugin-(flow|no-implicit-side-effects|consistent-subscribe)$/.test(pkgData.name)
		).filter(
			pkgData => /^https?:\/\/github.com\//.test(pkgData.links.repository)
		).map(
			pkgData => get(
				pkgData.links.repository + '/tree/HEAD/docs/rules',
				'.files a[href$=".md"]'
			).then(
				rules => rules.map(
					rule => rule.replace(/\/blob\/[^/]+/, '/blob/HEAD') + '#readme'
				)
			).then(
				rules => rules.length ? rules : got(
					pkgData.links.repository.replace(/\/\/github.com\//i, '//raw.githubusercontent.com/') + '/master/README.md'
				).then(
					res => {
						if (!/^(#+)\s+Rules$/im.test(res.body)) {
							return;
						}
						const pluginName = pkgData.name.slice(14);
						let baseUrl = EslintPluginDocBaseUrl[pluginName];
						if (!baseUrl) {
							baseUrl = rule => npmUrlPrefix + pluginName + '#' + rule;
						}
						const level = RegExp.$1;
						let md = RegExp.rightContext;
						md = md.slice(0, md.indexOf(RegExp('^' + level + '\\s+')));
						md = md.match(RegExp('^#' + level + '+\\s+.+?$', 'gm')) || [];
						md = md.map(
							rule => rule.replace(/<(\w+)>(.+?)(<\/\1>)/, '$2').replace(/<(\w+)>(.+?)(<\/\1>)/, '$2').replace(/^#+\s+(.+?)\s*$/, '$1').replace(/^.*\//, '')
						).filter(
							rule => /^[a-z]+(?:-[a-z]+)*$/.test(rule)
						).map(baseUrl);
						return md;
					},
					() => null
				)
			)
		)
	),

	// JSCS
	fs.readdir('node_modules/jscs/lib/rules').then(files => (
		files.filter(file => /\.js$/.test(file)).map(rule => (
			rule.replace(/\.\w+$/, '').replace(/-[\w]/g, char => char[1].toUpperCase())
		)).map(rule => (
			`http://jscs.info/rule/${rule}`
		))
	)),

	// CSSLint
	require('csslint').CSSLint.getRules().map(rule => rule.url).filter(Boolean),

	// TSLint
	fs.readdir('node_modules/tslint/lib/rules').then(files => (
		files.filter(file => /Rule\.js$/.test(file)).map(rule => (
			rule.replace(/Rule\.\w+$/, '').replace(/[A-Z]/g, char => '-' + char.toLowerCase())
		)).map(rule => (
			`https://palantir.github.io/tslint/rules/${rule}/`
		))
	)),

	// stylelint
	require('stylelint/lib/rules').map(rule => (
		`https://stylelint.io/user-guide/rules/${rule}/`
	)),

	// HTMLHint
	Object.keys(require('htmlhint').HTMLHint.rules).map(rule => (
		`https://github.com/yaniswang/HTMLHint/wiki/${rule}`
	)),
]).then(
	urls => awaitArray(urls)
).then(
	urls => expandArray(urls, [])
).then(
	urls => Promise.all([
		updateFile(
			'../lib/shorturl.json',
			urls,
			shortUrl
		),
		updateFile(
			'../lib/shorturl_cn.json',
			urls,
			shortUrlCn
		),
	])
).then(hasChange => {
	if (hasChange.some(Boolean)) {
		const spawnSync = require('child_process').spawnSync;
		const git = (args) => {
			return spawnSync(
				'git',
				args,
				{
					stdio: 'inherit',
				}
			);
		};
		git([
			'--no-pager',
			'diff',
			'--',
			'lib/*.json',
		]);
		git([
			'add',
			'lib/*.json',
		]);
		process.exitCode = 127;
	}
});

process.on('unhandledRejection', error => {
	console.error(error);
	process.exit(1);
});
