'use strict';

const messages = new Map([
	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/alt-require.js
	['An alt attribute must be present on <img> elements.', '<img>元素必須要有alt屬性。'],
	[/^The alt attribute of (.+?) must have a value\.$/, '$1必須要有alt屬性。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-lowercase.js
	[/^The attribute name of (.+?) must be in lowercase\.$/, '屬性名$1必須小寫。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-no-duplication.js
	[/^Duplicate of attribute name (.+?) was found\.$/, '發現重複的屬性$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-unsafe-chars.js
	[/^The value of attribute (.+?) cannot contain an unsafe char (.+?)\.$/, '屬性$1的值中不能包含不安全字元$2。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-value-double-quotes.js
	[/^The value of attribute (.+?) must be in double quotes\.$/, '屬性$1的值必須放在雙引號中。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-value-not-empty.js
	[/^The attribute (.+?) must have a value\.$/, '屬性$1必須有值。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/doctype-first.js
	['Doctype must be declared first.', 'doctype必須首先聲明。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/doctype-html5.js
	['Invalid doctype. Use: "<!DOCTYPE html>"', '無效的doctype，請使用：“<!DOCTYPE html>”。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/href-abs-or-rel.js
	[/^The value of the href attribute (.+?) must be absolute\.$/, '屬性href的值$1必須是絕對路徑。'],
	[/^The value of the href attribute (.+?) must be relative\.$/, '屬性href的值$1必須是相對路徑。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/id-class-ad-disabled.js
	[/^The value of attribute (.+?) cannot use the ad keyword\.$/, '屬性$1的值不能包含關鍵字“ad”。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/id-class-value.js
	['The id and class attribute values must be in lowercase and split by an underscore.', 'id和class屬性的值必須小寫，並由下劃線分割。'],
	['The id and class attribute values must be in lowercase and split by a dash.', 'id和class屬性的值必須小寫，並由連字元分割。'],
	['The id and class attribute values must meet the camelCase style.', 'id和class屬性的值必須符合camelCase風格。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/id-unique.js
	[/^The id value (.+?) must be unique\.$/, 'id的值$1必須是唯一的。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/inline-script-disabled.js
	[/^Inline script (.+?) cannot be used\.$/, '不能使用行內指令碼$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/inline-style-disabled.js
	[/^Inline style (.+?) cannot be used\.$/, '不能使用行內樣式$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/space-tab-mixed-disabled.js
	[/^Please use space for indentation and keep (.+?) length\.$/, '請使用$1個空格作為縮進。'],
	['Please use space for indentation.', '請使用空格作為縮進。'],
	['Please use tab for indentation.', '請使用tab作為縮進。'],
	['Do not mix tabs and spaces for indentation.', '不要混用tab和空格作為縮進。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/spec-char-escape.js
	[/^Special characters must be escaped ?: (.+?)\.$/, '必須轉義特殊字元：$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/src-not-empty.js
	[/^The attribute (.+?) of the tag (.+?) must have a value\.$/, '$1標籤的屬性$2必須有值。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/style-disabled.js
	['The <style> tag cannot be used.', '不能使用<style>標籤。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/tag-pair.js
	[/^Tag must be paired, missing: (.+?), (?:start|open) tag match failed (.+?) on line (.+?)\.$/, '標籤必須匹配，缺失：$1，在第$3行匹配開始標籤$2失敗'],
	[/^Tag must be paired, no start tag: (.+?)$/, '標籤必須匹配，沒有開始標籤：$1'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/tag-self-close.js
	[/^The empty tag ?: (.+?) must be self closed\.$/, '空標籤：$1必須自關閉'],

	[/^The html element name of (.+?) must be in lowercase\.$/, 'html元素$1的名字必須小寫。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/title-require.js
	['<title></title> must not be empty.', '<title></title> 必須不為空。'],
	['<title> must be present in <head> tag.', '<title>必須存在於<head>標籤。'],
]);
module.exports = function (error) {
	for (const [en, han] of messages) {
		if (en instanceof RegExp) {
			if (en.test(error.message)) {
				error.message = error.message.replace(en, han);
				return error;
			}
		} else if (error.message === en) {
			error.message = han;
			return error;
		}
	}
	return error;
};
