// Initial Template for C++
#include <bits/stdc++.h>
// using namespace std;

// User function template for C++
class Solution {
public:
    std::vector<int> factorial(int N) {
        std::vector<int> result;
        result.push_back(1);
        for (int x = 2; x <= N; x++) {
            int carry = 0;
            for (int i = 0; i < (int)result.size(); i++) {
                int prod = result[i] * x + carry;
                result[i] = prod % 10;
                carry = prod / 10;
            }
            while (carry != 0) {
                result.push_back(carry % 10);
                carry /= 10;
            }
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    std::vector<int> plusOne(std::vector<int>& digits) {
        // 傳入一個以一維vector(字串)儲存的非負整數 求該數字+1的值
        
        // 從最右邊 最後一位(LSB)開始檢查
        for (int i = digits.size() - 1; i >= 0; i--) {
            if (digits[i] == 9) // 如果最後一位是9 代表有進位問題
                digits[i] = 0;  // 給0(進位) 再繼續檢查前一位
            else if (digits[i] != 9) { // 如果最後一位不是9 則沒有進位問題
                digits[i] += 1; // 直接+1後 回傳結果
                return digits;
            }
        } 
        // check if there's any leading zero
        if (digits.front() == 0) {
            digits.insert(digits.begin(), 1);
        }
        return digits;
    }
};

// { Driver Code Starts.
int main() {
    int t;
    std::cin >> t;
    std::cout << "number of instances: " << t << "\n";
    while (t--) {
        int N;
        std::cin >> N;
        Solution ob;
        std::vector<int> result = ob.factorial(N);
        std::cout << "factorial of " << N << ": ";
        for (int i = 0; i < (int)result.size(); ++i){
            std::cout<< result[i];
        }
        std::cout << "\nnumber of digits: " << result.size() << "\n";
    }
    std::cout << std::endl;
    // int N = 200;
    // Solution ob;
    // std::vector<int> result = ob.factorial(N);
    // std::cout << "\nfactorial of " << N << ": ";
    // for (int i = 0; i < (int)result.size(); ++i){
    //     std::cout<< result[i];
    // }
    // std::cout << std::endl;
    // std::cout << "number of digits: " << result.size() << "\n\n";
    return 0;
}  // } Driver Code Ends
