#include <vector>
class IntStack
{
    public:
        /// Returns true if there are no elements on the stack
        bool empty() const{
            if (m_size==0){
                return true;
            }
            return false;
        }
        
        /// Returns the number of elements on the stack
        int size() const{
            return m_size;
        }

        /// Return the top element of the stack, i.e. the last element that was pushed
        int top() const{
            return m_vector[m_size-1];
        }

        /// Adds a new element to the top of the stack
        void push(const int& new_value){
            m_vector[m_size]=new_value;
            m_size++;
        }

        /// Remove the topmost element of the stack
        void pop(){
            m_vector.pop_back();     
            m_size--;
        }

    private:
        std::vector<int> m_vector = std::vector<int>(300000000); 
        int m_size = 0;
};

